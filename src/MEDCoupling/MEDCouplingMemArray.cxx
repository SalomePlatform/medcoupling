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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingMemArray.txx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "GenMathFormulae.hxx"
#include "InterpKernelExprParser.hxx"

#include <set>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <functional>

typedef double (*MYFUNCPTR)(double);

using namespace ParaMEDMEM;

template<int SPACEDIM>
void DataArrayDouble::findCommonTuplesAlg(const double *bbox, int nbNodes, int limitNodeId, double prec, DataArrayInt *c, DataArrayInt *cI) const
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
                  cI->pushBackSilent(cI->back()+(int)commonNodes.size()+1);
                  c->pushBackSilent(i);
                  c->insertAtTheEnd(commonNodes.begin(),commonNodes.end());
                }
            }
        }
    }
}

template<int SPACEDIM>
void DataArrayDouble::FindTupleIdsNearTuplesAlg(const BBTree<SPACEDIM,int>& myTree, const double *pos, int nbOfTuples, double eps,
                                                DataArrayInt *c, DataArrayInt *cI)
{
  for(int i=0;i<nbOfTuples;i++)
    {
      std::vector<int> intersectingElems;
      myTree.getElementsAroundPoint(pos+i*SPACEDIM,intersectingElems);
      std::vector<int> commonNodes;
      for(std::vector<int>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
        commonNodes.push_back(*it);
      cI->pushBackSilent(cI->back()+(int)commonNodes.size());
      c->insertAtTheEnd(commonNodes.begin(),commonNodes.end());
    }
}

template<int SPACEDIM>
void DataArrayDouble::FindClosestTupleIdAlg(const BBTreePts<SPACEDIM,int>& myTree, double dist, const double *pos, int nbOfTuples, const double *thisPt, int thisNbOfTuples, int *res)
{
  double distOpt(dist);
  const double *p(pos);
  int *r(res);
  for(int i=0;i<nbOfTuples;i++,p+=SPACEDIM,r++)
    {
      while(true)
        {
          int elem=-1;
          double ret=myTree.getElementsAroundPoint2(p,distOpt,elem);
          if(ret!=std::numeric_limits<double>::max())
            {
              distOpt=std::max(ret,1e-4);
              *r=elem;
              break;
            }
          else
            { distOpt=2*distOpt; continue; }
        }
    }
}

std::size_t DataArray::getHeapMemorySize() const
{
  std::size_t sz1=_name.capacity();
  std::size_t sz2=_info_on_compo.capacity();
  std::size_t sz3=0;
  for(std::vector<std::string>::const_iterator it=_info_on_compo.begin();it!=_info_on_compo.end();it++)
    sz3+=(*it).capacity();
  return sz1+sz2+sz3;
}

/*!
 * Sets the attribute \a _name of \a this array.
 * See \ref MEDCouplingArrayBasicsName "DataArrays infos" for more information.
 *  \param [in] name - new array name
 */
void DataArray::setName(const char *name)
{
  _name=name;
}

/*!
 * Copies textual data from an \a other DataArray. The copied data are
 * - the name attribute,
 * - the information of components.
 *
 * For more information on these data see \ref MEDCouplingArrayBasicsName "DataArrays infos".
 *
 *  \param [in] other - another instance of DataArray to copy the textual data from.
 *  \throw If number of components of \a this array differs from that of the \a other.
 */
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

/*!
 * Compares textual information of \a this DataArray with that of an \a other one.
 * The compared data are
 * - the name attribute,
 * - the information of components.
 *
 * For more information on these data see \ref MEDCouplingArrayBasicsName "DataArrays infos".
 *  \param [in] other - another instance of DataArray to compare the textual data of.
 *  \return bool - \a true if the textual information is same, \a false else.
 */
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

std::string DataArray::cppRepr(const char *varName) const throw(INTERP_KERNEL::Exception)
{
  std::ostringstream ret;
  reprCppStream(varName,ret);
  return ret.str();
}

/*!
 * Sets information on all components. To know more on format of this information
 * see \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - a vector of strings.
 *  \throw If size of \a info differs from the number of components of \a this.
 */
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

/*!
 * Returns information on a component specified by an index.
 * To know more on format of this information
 * see \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \return std::string - a string containing the information on \a i-th component.
 *  \throw If \a i is not a valid component index.
 */
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
 * Returns the var part of the full information of the \a i-th component.
 * For example, if \c getInfoOnComponent(0) returns "SIGXY [N/m^2]", then
 * \c getVarOnComponent(0) returns "SIGXY".
 * If a unit part of information is not detected by presence of
 * two square brackets, then the full information is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \return std::string - a string containing the var information, or the full info.
 *  \throw If \a i is not a valid component index.
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
 * Returns the unit part of the full information of the \a i-th component.
 * For example, if \c getInfoOnComponent(0) returns "SIGXY [ N/m^2]", then
 * \c getUnitOnComponent(0) returns " N/m^2".
 * If a unit part of information is not detected by presence of
 * two square brackets, then an empty string is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \return std::string - a string containing the unit information, if any, or "".
 *  \throw If \a i is not a valid component index.
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

/*!
 * Returns the var part of the full component information.
 * For example, if \a info == "SIGXY [N/m^2]", then this method returns "SIGXY".
 * If a unit part of information is not detected by presence of
 * two square brackets, then the whole \a info is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - the full component information.
 *  \return std::string - a string containing only var information, or the \a info.
 */
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

/*!
 * Returns the unit part of the full component information.
 * For example, if \a info == "SIGXY [ N/m^2]", then this method returns " N/m^2".
 * If a unit part of information is not detected by presence of
 * two square brackets, then an empty string is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - the full component information.
 *  \return std::string - a string containing only unit information, if any, or "".
 */
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

/*!
 * Sets information on a component specified by an index.
 * To know more on format of this information
 * see \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \param [in] info - the string containing the information.
 *  \throw If \a i is not a valid component index.
 *  \warning Don't pass NULL as \a info!
 */
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
          std::ostringstream oss; oss << "DataArray::CheckValueInRangeEx : " << msg  << " ! Expected start " << start << " of input range, in [0," << value << "[ !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(end<0 || end>value)
    {
      std::ostringstream oss; oss << "DataArray::CheckValueInRangeEx : " << msg  << " ! Expected end " << end << " of input range, in [0," << value << "] !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::CheckClosingParInRange(int ref, int value, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(value<0 || value>ref)
    {
      std::ostringstream oss; oss << "DataArray::CheckClosingParInRange : " << msg  << " ! Expected input range in [0," << ref << "] having closing open parenthesis " << value << " !";
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
  if(end==begin)
    return 0;
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

/*!
 * Returns a new instance of DataArrayDouble. The caller is to delete this array
 * using decrRef() as it is no more needed. 
 */
DataArrayDouble *DataArrayDouble::New()
{
  return new DataArrayDouble;
}

/*!
 * Checks if raw data is allocated. Read more on the raw data
 * in \ref MEDCouplingArrayBasicsTuplesAndCompo "DataArrays infos" for more information.
 *  \return bool - \a true if the raw data is allocated, \a false else.
 */
bool DataArrayDouble::isAllocated() const
{
  return getConstPointer()!=0;
}

/*!
 * Checks if raw data is allocated and throws an exception if it is not the case.
 *  \throw If the raw data is not allocated.
 */
void DataArrayDouble::checkAllocated() const throw(INTERP_KERNEL::Exception)
{
  if(!isAllocated())
    throw INTERP_KERNEL::Exception("DataArrayDouble::checkAllocated : Array is defined but not allocated ! Call alloc or setValues method first !");
}

std::size_t DataArrayDouble::getHeapMemorySize() const
{
  std::size_t sz=(std::size_t)_mem.getNbOfElemAllocated();
  sz*=sizeof(double);
  return DataArray::getHeapMemorySize()+sz;
}

/*!
 * Sets information on all components. This method can change number of components
 * at certain conditions; if the conditions are not respected, an exception is thrown.
 * The number of components can be changed provided that \a this is not allocated.
 *
 * To know more on format of the component information see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - a vector of component infos.
 *  \throw If \a this->getNumberOfComponents() != \a info.size() && \a this->isAllocated()
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
 * Returns the only one value in \a this, if and only if number of elements
 * (nb of tuples * nb of components) is equal to 1, and that \a this is allocated.
 *  \return double - the sole value stored in \a this array.
 *  \throw If at least one of conditions stated above is not fulfilled.
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
 * Checks the number of tuples.
 *  \return bool - \a true if getNumberOfTuples() == 0, \a false else.
 *  \throw If \a this is not allocated.
 */
bool DataArrayDouble::empty() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  return getNumberOfTuples()==0;
}

/*!
 * Returns a full copy of \a this. For more info on copying data arrays see
 * \ref MEDCouplingArrayBasicsCopyDeep.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble.
 */
DataArrayDouble *DataArrayDouble::deepCpy() const
{
  return new DataArrayDouble(*this);
}

/*!
 * Returns either a \a deep or \a shallow copy of this array. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeep and \ref MEDCouplingArrayBasicsCopyShallow.
 *  \param [in] dCpy - if \a true, a deep copy is returned, else, a shallow one.
 *  \return DataArrayDouble * - either a new instance of DataArrayDouble (if \a dCpy
 *          == \a true) or \a this instance (if \a dCpy == \a false).
 */
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

/*!
 * Copies all the data from another DataArrayDouble. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeepAssign.
 *  \param [in] other - another instance of DataArrayDouble to copy data from.
 *  \throw If the \a other is not allocated.
 */
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

/*!
 * This method reserve nbOfElems elements in memory ( nbOfElems*8 bytes ) \b without impacting the number of tuples in \a this.
 * If \a this has already been allocated, this method checks that \a this has only one component. If not an INTERP_KERNEL::Exception will be thrown.
 * If \a this has not already been allocated, number of components is set to one.
 * This method allows to reduce number of reallocations on invokation of DataArrayDouble::pushBackSilent and DataArrayDouble::pushBackValsSilent on \a this.
 * 
 * \sa DataArrayDouble::pack, DataArrayDouble::pushBackSilent, DataArrayDouble::pushBackValsSilent
 */
void DataArrayDouble::reserve(int nbOfElems) throw(INTERP_KERNEL::Exception)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    {
      _mem.reserve(nbOfElems);
    }
  else if(nbCompo==0)
    {
      _mem.reserve(nbOfElems);
      _info_on_compo.resize(1);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayDouble::reserve : not available for DataArrayDouble with number of components different than 1 !");
}

/*!
 * This method adds at the end of \a this the single value \a val. This method do \b not update its time label to avoid useless incrementation
 * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
 *
 * \param [in] val the value to be added in \a this
 * \throw If \a this has already been allocated with number of components different from one.
 * \sa DataArrayDouble::pushBackValsSilent
 */
void DataArrayDouble::pushBackSilent(double val) throw(INTERP_KERNEL::Exception)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    _mem.pushBack(val);
  else if(nbCompo==0)
    {
      _info_on_compo.resize(1);
      _mem.pushBack(val);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayDouble::pushBackSilent : not available for DataArrayDouble with number of components different than 1 !");
}

/*!
 * This method adds at the end of \a this a serie of values [\c valsBg,\c valsEnd). This method do \b not update its time label to avoid useless incrementation
 * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
 *
 *  \param [in] valsBg - an array of values to push at the end of \this.
 *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
 *              the last value of \a valsBg is \a valsEnd[ -1 ].
 * \throw If \a this has already been allocated with number of components different from one.
 * \sa DataArrayDouble::pushBackSilent
 */
void DataArrayDouble::pushBackValsSilent(const double *valsBg, const double *valsEnd) throw(INTERP_KERNEL::Exception)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    _mem.insertAtTheEnd(valsBg,valsEnd);
  else if(nbCompo==0)
    {
      _info_on_compo.resize(1);
      _mem.insertAtTheEnd(valsBg,valsEnd);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayDouble::pushBackValsSilent : not available for DataArrayDouble with number of components different than 1 !");
}

/*!
 * This method returns silently ( without updating time label in \a this ) the last value, if any and suppress it.
 * \throw If \a this is already empty.
 * \throw If \a this has number of components different from one.
 */
double DataArrayDouble::popBackSilent() throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()==1)
    return _mem.popBack();
  else
    throw INTERP_KERNEL::Exception("DataArrayDouble::popBackSilent : not available for DataArrayDouble with number of components different than 1 !");
}

/*!
 * This method \b do \b not modify content of \a this. It only modify its memory footprint if the allocated memory is to high regarding real data to store.
 *
 * \sa DataArrayDouble::getHeapMemorySize, DataArrayDouble::reserve
 */
void DataArrayDouble::pack() const throw(INTERP_KERNEL::Exception)
{
  _mem.pack();
}

/*!
 * Allocates the raw data in memory. If exactly same memory as needed already
 * allocated, it is not re-allocated.
 *  \param [in] nbOfTuple - number of tuples of data to allocate.
 *  \param [in] nbOfCompo - number of components of data to allocate.
 *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
 */
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

/*!
 * Allocates the raw data in memory. If the memory was already allocated, then it is
 * freed and re-allocated. See an example of this method use
 * \ref MEDCouplingArraySteps1WC "here".
 *  \param [in] nbOfTuple - number of tuples of data to allocate.
 *  \param [in] nbOfCompo - number of components of data to allocate.
 *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
 */
void DataArrayDouble::alloc(int nbOfTuple, int nbOfCompo) throw(INTERP_KERNEL::Exception)
{
  if(nbOfTuple<0 || nbOfCompo<0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::alloc : request for negative length of data !");
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*nbOfTuple);
  declareAsNew();
}

/*!
 * Assign zero to all values in \a this array. To know more on filling arrays see
 * \ref MEDCouplingArrayFill.
 * \throw If \a this is not allocated.
 */
void DataArrayDouble::fillWithZero() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(0.);
  declareAsNew();
}

/*!
 * Assign \a val to all values in \a this array. To know more on filling arrays see
 * \ref MEDCouplingArrayFill.
 *  \param [in] val - the value to fill with.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::fillWithValue(double val) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(val);
  declareAsNew();
}

/*!
 * Set all values in \a this array so that the i-th element equals to \a init + i
 * (i starts from zero). To know more on filling arrays see \ref MEDCouplingArrayFill.
 *  \param [in] init - value to assign to the first element of array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Checks if all values in \a this array are equal to \a val at precision \a eps.
 *  \param [in] val - value to check equality of array values to.
 *  \param [in] eps - precision to check the equality.
 *  \return bool - \a true if all values are in range (_val_ - _eps_; _val_ + _eps_),
 *                 \a false else.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Sorts values of the array.
 *  \param [in] asc - \a true means ascending order, \a false, descending.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
void DataArrayDouble::sort(bool asc) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::sort : only supported with 'this' array with ONE component !");
  _mem.sort(asc);
}

/*!
 * Reverse the array values.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::reverse() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::reverse : only supported with 'this' array with ONE component !");
  _mem.reverse();
}

/*!
 * Checks that \a this array is consistently **increasing** or **decreasing** in value,
 * with at least absolute difference value of |\a eps| at each step.
 * If not an exception is thrown.
 *  \param [in] increasing - if \a true, the array values should be increasing.
 *  \param [in] eps - minimal absolute difference between the neighbor values at which 
 *                    the values are considered different.
 *  \throw If sequence of values is not strictly monotonic in agreement with \a
 *         increasing arg.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::checkMonotonic(bool increasing, double eps) const throw(INTERP_KERNEL::Exception)
{
  if(!isMonotonic(increasing,eps))
    {
      if (increasing)
        throw INTERP_KERNEL::Exception("DataArrayDouble::checkMonotonic : 'this' is not INCREASING monotonic !");
      else
        throw INTERP_KERNEL::Exception("DataArrayDouble::checkMonotonic : 'this' is not DECREASING monotonic !");
    }
}

/*!
 * Checks that \a this array is consistently **increasing** or **decreasing** in value,
 * with at least absolute difference value of |\a eps| at each step.
 *  \param [in] increasing - if \a true, array values should be increasing.
 *  \param [in] eps - minimal absolute difference between the neighbor values at which 
 *                    the values are considered different.
 *  \return bool - \a true if values change in accordance with \a increasing arg.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
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
  if(increasing)
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

/*!
 * Returns a textual and human readable representation of \a this instance of
 * DataArrayDouble. This text is shown when a DataArrayDouble is printed in Python.
 *  \return std::string - text describing \a this DataArrayDouble.
 */
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
  ofs.precision(17);
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
  stream.precision(17);
  _mem.repr(getNumberOfComponents(),stream);
}

void DataArrayDouble::reprZipWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  stream.precision(17);
  _mem.reprZip(getNumberOfComponents(),stream);
}

void DataArrayDouble::reprCppStream(const char *varName, std::ostream& stream) const
{
  int nbTuples=getNumberOfTuples(),nbComp=getNumberOfComponents();
  const double *data=getConstPointer();
  stream.precision(17);
  stream << "DataArrayDouble *" << varName << "=DataArrayDouble::New();" << std::endl;
  if(nbTuples*nbComp>=1)
    {
      stream << "const double " << varName << "Data[" << nbTuples*nbComp << "]={";
      std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<double>(stream,","));
      stream << data[nbTuples*nbComp-1] << "};" << std::endl;
      stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
    }
  else
    stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
  stream << varName << "->setName(\"" << getName() << "\");" << std::endl;
}

bool DataArrayDouble::isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,prec,reason);
}

/*!
 * Checks if \a this and another DataArrayDouble are fully equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayDouble to compare with \a this one.
 *  \param [in] prec - precision value to compare numeric data of the arrays.
 *  \return bool - \a true if the two arrays are equal, \a false else.
 */
bool DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,prec,tmp);
}

/*!
 * Checks if values of \a this and another DataArrayDouble are equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayDouble to compare with \a this one.
 *  \param [in] prec - precision value to compare numeric data of the arrays.
 *  \return bool - \a true if the values of two arrays are equal, \a false else.
 */
bool DataArrayDouble::isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,prec,tmp);
}

/*!
 * Changes number of tuples in the array. If the new number of tuples is smaller
 * than the current number the array is truncated, otherwise the array is extended.
 *  \param [in] nbOfTuples - new number of tuples. 
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.reAlloc(getNumberOfComponents()*nbOfTuples);
  declareAsNew();
}

/*!
 * Creates a new DataArrayInt and assigns all (textual and numerical) data of \a this
 * array to the new one.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 */
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

/*!
 * Returns a new DataArrayDouble holding the same values as \a this array but differently
 * arranged in memory. If \a this array holds 2 components of 3 values:
 * \f$ x_0,x_1,x_2,y_0,y_1,y_2 \f$, then the result array holds these values arranged
 * as follows: \f$ x_0,y_0,x_1,y_1,x_2,y_2 \f$.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \warning Do not confuse this method with transpose()!
 */
DataArrayDouble *DataArrayDouble::fromNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromNoInterlace : Not defined array !");
  double *tab=_mem.fromNoInterlace(getNumberOfComponents());
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

/*!
 * Returns a new DataArrayDouble holding the same values as \a this array but differently
 * arranged in memory. If \a this array holds 2 components of 3 values:
 * \f$ x_0,y_0,x_1,y_1,x_2,y_2 \f$, then the result array holds these values arranged
 * as follows: \f$ x_0,x_1,x_2,y_0,y_1,y_2 \f$.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \warning Do not confuse this method with transpose()!
 */
DataArrayDouble *DataArrayDouble::toNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayDouble::toNoInterlace : Not defined array !");
  double *tab=_mem.toNoInterlace(getNumberOfComponents());
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

/*!
 * Permutes values of \a this array as required by \a old2New array. The values are
 * permuted so that \c new[ \a old2New[ i ]] = \c old[ i ]. Number of tuples remains
 * the same as in \this one.
 * If a permutation reduction is needed, substr() or selectByTupleId() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *     giving a new position for i-th old value.
 */
void DataArrayDouble::renumberInPlace(const int *old2New)
{
  checkAllocated();
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
 * Permutes values of \a this array as required by \a new2Old array. The values are
 * permuted so that \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of tuples remains
 * the same as in \this one.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
 *     giving a previous position of i-th new value.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
void DataArrayDouble::renumberInPlaceR(const int *new2Old)
{
  checkAllocated();
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
 * Returns a copy of \a this array with values permuted as required by \a old2New array.
 * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ].
 * Number of tuples in the result array remains the same as in \this one.
 * If a permutation reduction is needed, renumberAndReduce() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *          giving a new position for i-th old value.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 */
DataArrayDouble *DataArrayDouble::renumber(const int *old2New) const
{
  checkAllocated();
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
 * Returns a copy of \a this array with values permuted as required by \a new2Old array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of
 * tuples in the result array remains the same as in \this one.
 * If a permutation reduction is needed, substr() or selectByTupleId() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
 *     giving a previous position of i-th new value.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayDouble *DataArrayDouble::renumberR(const int *new2Old) const
{
  checkAllocated();
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
 * Returns a shorten and permuted copy of \a this array. The new DataArrayDouble is
 * of size \a newNbOfTuple and it's values are permuted as required by \a old2New array.
 * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ] for all
 * \a old2New[ i ] >= 0. In other words every i-th tuple in \a this array, for which 
 * \a old2New[ i ] is negative, is missing from the result array.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *     giving a new position for i-th old tuple and giving negative position for
 *     for i-th old tuple that should be omitted.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayDouble *DataArrayDouble::renumberAndReduce(const int *old2New, int newNbOfTuple) const
{
  checkAllocated();
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
 * Returns a shorten and permuted copy of \a this array. The new DataArrayDouble is
 * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
 * \a new2OldBg array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
 * This method is equivalent to renumberAndReduce() except that convention in input is
 * \c new2old and \b not \c old2new.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
 *              tuple index in \a this array to fill the i-th tuple in the new array.
 *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
 *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
 *              \a new2OldBg <= \a pi < \a new2OldEnd.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayDouble *DataArrayDouble::selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  double *pt=ret->getPointer();
  const double *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten and permuted copy of \a this array. The new DataArrayDouble is
 * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
 * \a new2OldBg array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
 * This method is equivalent to renumberAndReduce() except that convention in input is
 * \c new2old and \b not \c old2new.
 * This method is equivalent to selectByTupleId() except that it prevents coping data
 * from behind the end of \a this array.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
 *              tuple index in \a this array to fill the i-th tuple in the new array.
 *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
 *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
 *              \a new2OldBg <= \a pi < \a new2OldEnd.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a new2OldEnd - \a new2OldBg > \a this->getNumberOfTuples().
 */
DataArrayDouble *DataArrayDouble::selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
  return ret.retn();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayDouble contains every
 * (\a bg + \c i * \a step)-th tuple of \a this array located before the \a end2-th
 * tuple. Indices of the selected tuples are the same as ones returned by the Python
 * command \c range( \a bg, \a end2, \a step ).
 * This method is equivalent to selectByTupleIdSafe() except that the input array is
 * not constructed explicitly.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] bg - index of the first tuple to copy from \a this array.
 *  \param [in] end2 - index of the tuple before which the tuples to copy are located.
 *  \param [in] step - index increment to get index of the next tuple to copy.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If (\a end2 < \a bg) or (\a step <= 0).
 *  \sa DataArrayDouble::substr.
 */
DataArrayDouble *DataArrayDouble::selectByTupleId2(int bg, int end2, int step) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  int newNbOfTuples=GetNumberOfItemGivenBES(bg,end2,step,"DataArrayDouble::selectByTupleId2 : ");
  ret->alloc(newNbOfTuples,nbComp);
  double *pt=ret->getPointer();
  const double *srcPt=getConstPointer()+bg*nbComp;
  for(int i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
    std::copy(srcPt,srcPt+nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayDouble contains ranges
 * of tuples specified by \a ranges parameter.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] ranges - std::vector of std::pair's each of which defines a range
 *              of tuples in [\c begin,\c end) format.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a end < \a begin.
 *  \throw If \a end > \a this->getNumberOfTuples().
 *  \throw If \a this is not allocated.
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
  return ret.retn();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayDouble contains all
 * tuples starting from the \a tupleIdBg-th tuple and including all tuples located before
 * the \a tupleIdEnd-th one. This methods has a similar behavior as std::string::substr().
 * This method is a specialization of selectByTupleId2().
 *  \param [in] tupleIdBg - index of the first tuple to copy from \a this array.
 *  \param [in] tupleIdEnd - index of the tuple before which the tuples to copy are located.
 *          If \a tupleIdEnd == -1, all the tuples till the end of \a this array are copied.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a tupleIdBg < 0.
 *  \throw If \a tupleIdBg > \a this->getNumberOfTuples().
    \throw If \a tupleIdEnd != -1 && \a tupleIdEnd < \a this->getNumberOfTuples().
 *  \sa DataArrayDouble::selectByTupleId2
 */
DataArrayDouble *DataArrayDouble::substr(int tupleIdBg, int tupleIdEnd) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Returns a shorten or extended copy of \a this array. If \a newNbOfComp is less
 * than \a this->getNumberOfComponents() then the result array is shorten as each tuple
 * is truncated to have \a newNbOfComp components, keeping first components. If \a
 * newNbOfComp is more than \a this->getNumberOfComponents() then the result array is
 * expanded as each tuple is populated with \a dftValue to have \a newNbOfComp
 * components.  
 *  \param [in] newNbOfComp - number of components for the new array to have.
 *  \param [in] dftValue - value assigned to new values added to the new array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
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
 * Changes the number of components within \a this array so that its raw data **does
 * not** change, instead splitting this data into tuples changes.
 *  \param [in] newNbOfComp - number of components for \a this array to have.
 *  \throw If \a this is not allocated
 *  \throw If getNbOfElems() % \a newNbOfCompo != 0.
 *  \warning This method erases all (name and unit) component info set before!
 */
void DataArrayDouble::rearrange(int newNbOfCompo) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfElems=getNbOfElems();
  if(nbOfElems%newNbOfCompo!=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::rearrange : nbOfElems%newNbOfCompo!=0 !");
  _info_on_compo.clear();
  _info_on_compo.resize(newNbOfCompo);
  declareAsNew();
}

/*!
 * Changes the number of components within \a this array to be equal to its number
 * of tuples, and inversely its number of tuples to become equal to its number of 
 * components. So that its raw data **does not** change, instead splitting this
 * data into tuples changes.
 *  \throw If \a this is not allocated.
 *  \warning This method erases all (name and unit) component info set before!
 *  \warning Do not confuse this method with fromNoInterlace() and toNoInterlace()!
 *  \sa rearrange()
 */
void DataArrayDouble::transpose() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  rearrange(nbOfTuples);
}

/*!
 * Returns a copy of \a this array composed of selected components.
 * The new DataArrayDouble has the same number of tuples but includes components
 * specified by \a compoIds parameter. So that getNbOfElems() of the result array
 * can be either less, same or more than \a this->getNbOfElems().
 *  \param [in] compoIds - sequence of zero based indices of components to include
 *              into the new array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If a component index (\a i) is not valid: 
 *         \a i < 0 || \a i >= \a this->getNumberOfComponents().
 *
 *  \ref cpp_mcdataarraydouble_keepselectedcomponents "Here is a Python example".
 */
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
  return ret.retn();
}

/*!
 * Appends components of another array to components of \a this one, tuple by tuple.
 * So that the number of tuples of \a this array remains the same and the number of 
 * components increases.
 *  \param [in] other - the DataArrayDouble to append to \a this one.
 *  \throw If \a this is not allocated.
 *  \throw If \a this and \a other arrays have different number of tuples.
 *
 *  \ref cpp_mcdataarraydouble_meldwith "Here is a C++ example".
 *
 *  \ref py_mcdataarraydouble_meldwith "Here is a Python example".
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
 * Searches for tuples coincident within \a prec tolerance. Each tuple is considered
 * as coordinates of a point in getNumberOfComponents()-dimensional space. The
 * distance is computed using norm2.
 *
 * Indices of coincident tuples are stored in output arrays.
 * A pair of arrays (\a comm, \a commIndex) is called "Surjective Format 2".
 *
 * This method is typically used by MEDCouplingPointSet::findCommonNodes() and
 * MEDCouplingUMesh::mergeNodes().
 *  \param [in] prec - minimal absolute distance between two tuples at which they are
 *              considered not coincident.
 *  \param [in] limitTupleId - limit tuple id. Tuples with id strictly lower than \a 
 *              limitTupleId are not considered.
 *  \param [out] comm - the array holding ids (== indices) of coincident tuples. 
 *               \a comm->getNumberOfComponents() == 1. 
 *               \a comm->getNumberOfTuples() == \a commIndex->back().
 *  \param [out] commIndex - the array dividing all indices stored in \a comm into
 *               groups of (indices of) coincident tuples. Its every value is a tuple
 *               index where a next group of tuples begins. For example the second
 *               group of tuples in \a comm is described by following range of indices:
 *               [ \a commIndex[1], \a commIndex[2] ). \a commIndex->getNumberOfTuples()-1
 *               gives the number of groups of coincident tuples.
 *  \throw If \a this is not allocated.
 *  \throw If \a this and \a other arrays have different number of tuples.
 *  \throw If the number of components is not in [1,2,3].
 *
 *  \ref cpp_mcdataarraydouble_findcommontuples "Here is a C++ example".
 *
 *  \ref py_mcdataarraydouble_findcommontuples  "Here is a Python example".
 *  \sa DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2().
 */
void DataArrayDouble::findCommonTuples(double prec, int limitTupleId, DataArrayInt *&comm, DataArrayInt *&commIndex) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if ((nbOfCompo<1) || (nbOfCompo>3)) //test before work
    throw INTERP_KERNEL::Exception("DataArrayDouble::findCommonTuples : Unexpected spacedim of coords. Must be 1, 2 or 3.");
  
  int nbOfTuples=getNumberOfTuples();
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> bbox=computeBBoxPerTuple(prec);
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c(DataArrayInt::New()),cI(DataArrayInt::New()); c->alloc(0,1); cI->pushBackSilent(0);
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
  comm=c.retn();
  commIndex=cI.retn();
}

/*!
 * 
 * \param [in] nbTimes specifies the nb of times each tuples in \a this will be duplicated contiguouly in returned DataArrayDouble instance.
 *             \a nbTimes  should be at least equal to 1.
 * \return a newly allocated DataArrayDouble having one component and number of tuples equal to \a nbTimes * \c this->getNumberOfTuples.
 * \throw if \a this is not allocated or if \a this has not number of components set to one or if \a nbTimes is lower than 1.
 */
DataArrayDouble *DataArrayDouble::duplicateEachTupleNTimes(int nbTimes) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::duplicateEachTupleNTimes : this should have only one component !");
  if(nbTimes<1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::duplicateEachTupleNTimes : nb times should be >= 1 !");
  int nbTuples=getNumberOfTuples();
  const double *inPtr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New(); ret->alloc(nbTimes*nbTuples,1);
  double *retPtr=ret->getPointer();
  for(int i=0;i<nbTuples;i++,inPtr++)
    {
      double val=*inPtr;
      for(int j=0;j<nbTimes;j++,retPtr++)
        *retPtr=val;
    }
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * This methods returns the minimal distance between the two set of points \a this and \a other.
 * So \a this and \a other have to have the same number of components. If not an INTERP_KERNEL::Exception will be thrown.
 * This method works only if number of components of \a this (equal to those of \a other) is in 1, 2 or 3.
 *
 * \param [out] thisTupleId the tuple id in \a this corresponding to the returned minimal distance
 * \param [out] otherTupleId the tuple id in \a other corresponding to the returned minimal distance
 * \return the minimal distance between the two set of points \a this and \a other.
 * \sa DataArrayDouble::findClosestTupleId
 */
double DataArrayDouble::minimalDistanceTo(const DataArrayDouble *other, int& thisTupleId, int& otherTupleId) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> part1=findClosestTupleId(other);
  int nbOfCompo(getNumberOfComponents());
  int otherNbTuples(other->getNumberOfTuples());
  const double *thisPt(begin()),*otherPt(other->begin());
  const int *part1Pt(part1->begin());
  double ret=std::numeric_limits<double>::max();
  for(int i=0;i<otherNbTuples;i++,part1Pt++,otherPt+=nbOfCompo)
    {
      double tmp(0.);
      for(int j=0;j<nbOfCompo;j++)
        tmp+=(otherPt[j]-thisPt[nbOfCompo*(*part1Pt)+j])*(otherPt[j]-thisPt[nbOfCompo*(*part1Pt)+j]);
      if(tmp<ret)
        { ret=tmp; thisTupleId=*part1Pt; otherTupleId=i; }
    }
  return sqrt(ret);
}

/*!
 * This methods returns for each tuple in \a other which tuple in \a this is the closest.
 * So \a this and \a other have to have the same number of components. If not an INTERP_KERNEL::Exception will be thrown.
 * This method works only if number of components of \a this (equal to those of \a other) is in 1, 2 or 3.
 *
 * \return a newly allocated (new object to be dealt by the caller) DataArrayInt having \c other->getNumberOfTuples() tuples and one components.
 * \sa DataArrayDouble::minimalDistanceTo
 */
DataArrayInt *DataArrayDouble::findClosestTupleId(const DataArrayDouble *other) const throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::findClosestTupleId : other instance is NULL !");
  checkAllocated(); other->checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo!=other->getNumberOfComponents())
    {
      std::ostringstream oss; oss << "DataArrayDouble::findClosestTupleId : number of components in this is " << nbOfCompo;
      oss << ", whereas number of components in other is " << other->getNumberOfComponents() << "! Should be equal !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfTuples=other->getNumberOfTuples();
  int thisNbOfTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(nbOfTuples,1);
  double bounds[6];
  getMinMaxPerComponent(bounds);
  switch(nbOfCompo)
    {
    case 3:
      {
        double xDelta(fabs(bounds[1]-bounds[0])),yDelta(fabs(bounds[3]-bounds[2])),zDelta(fabs(bounds[5]-bounds[4]));
        double delta=std::max(xDelta,yDelta); delta=std::max(delta,zDelta);
        double characSize=pow((delta*delta*delta)/((double)thisNbOfTuples),1./3.);
        BBTreePts<3,int> myTree(begin(),0,0,getNumberOfTuples(),characSize*1e-12);
        FindClosestTupleIdAlg<3>(myTree,3.*characSize*characSize,other->begin(),nbOfTuples,begin(),thisNbOfTuples,ret->getPointer());
        break;
      }
    case 2:
      {
        double xDelta(fabs(bounds[1]-bounds[0])),yDelta(fabs(bounds[3]-bounds[2]));
        double delta=std::max(xDelta,yDelta);
        double characSize=sqrt(delta/(double)thisNbOfTuples);
        BBTreePts<2,int> myTree(begin(),0,0,getNumberOfTuples(),characSize*1e-12);
        FindClosestTupleIdAlg<2>(myTree,2.*characSize*characSize,other->begin(),nbOfTuples,begin(),thisNbOfTuples,ret->getPointer());
        break;
      }
    case 1:
      {
        double characSize=fabs(bounds[1]-bounds[0])/thisNbOfTuples;
        BBTreePts<1,int> myTree(begin(),0,0,getNumberOfTuples(),characSize*1e-12);
        FindClosestTupleIdAlg<1>(myTree,1.*characSize*characSize,other->begin(),nbOfTuples,begin(),thisNbOfTuples,ret->getPointer());
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Unexpected spacedim of coords for findClosestTupleId. Must be 1, 2 or 3.");
    }
  return ret.retn();
}

/*!
 * Returns a copy of \a this array by excluding coincident tuples. Each tuple is
 * considered as coordinates of a point in getNumberOfComponents()-dimensional
 * space. The distance between tuples is computed using norm2. If several tuples are
 * not far each from other than \a prec, only one of them remains in the result
 * array. The order of tuples in the result array is same as in \a this one except
 * that coincident tuples are excluded.
 *  \param [in] prec - minimal absolute distance between two tuples at which they are
 *              considered not coincident.
 *  \param [in] limitTupleId - limit tuple id. Tuples with id strictly lower than \a 
 *              limiTupleId are not considered and thus not excluded.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If the number of components is not in [1,2,3].
 *
 *  \ref cpp_mcdataarraydouble_getdifferentvalues "Here is a Python example".
 */
DataArrayDouble *DataArrayDouble::getDifferentValues(double prec, int limitTupleId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayInt *c0=0,*cI0=0;
  findCommonTuples(prec,limitTupleId,c0,cI0);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c(c0),cI(cI0);
  int newNbOfTuples=-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(getNumberOfTuples(),c0->begin(),cI0->begin(),cI0->end(),newNbOfTuples);
  return renumberAndReduce(o2n->getConstPointer(),newNbOfTuples);
}

/*!
 * Copy all components in a specified order from another DataArrayDouble.
 * The specified components become the first ones in \a this array.
 * Both numerical and textual data is copied. The number of tuples in \a this and
 * the other array can be different.
 *  \param [in] a - the array to copy data from.
 *  \param [in] compoIds - sequence of zero based indices of components, data of which is
 *              to be copied.
 *  \throw If \a a is NULL.
 *  \throw If \a compoIds.size() != \a a->getNumberOfComponents().
 *  \throw If \a compoIds[i] < 0 or \a compoIds[i] > \a this->getNumberOfComponents().
 *
 *  \ref cpp_mcdataarraydouble_setselectedcomponents "Here is a Python example".
 */
void DataArrayDouble::setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setSelectedComponents : input DataArrayDouble is NULL !");
  checkAllocated();
  copyPartOfStringInfoFrom2(compoIds,*a);
  std::size_t partOfCompoSz=compoIds.size();
  int nbOfCompo=getNumberOfComponents();
  int nbOfTuples=std::min(getNumberOfTuples(),a->getNumberOfTuples());
  const double *ac=a->getConstPointer();
  double *nc=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(std::size_t j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}

/*!
 * Copy all values from another DataArrayDouble into specified tuples and components
 * of \a this array. Textual data is not copied.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign values to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - index of the first component of \a this array to assign values to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \param [in] strictCompoCompare - if \a true (by default), then \a a->getNumberOfComponents() 
 *              must be equal to the number of columns to assign to, else an
 *              exception is thrown; if \a false, then it is only required that \a
 *              a->getNbOfElems() equals to number of values to assign to (this condition
 *              must be respected even if \a strictCompoCompare is \a true). The number of 
 *              values to assign to is given by following Python expression:
 *              \a nbTargetValues = 
 *              \c len(\c range(\a bgTuples,\a endTuples,\a stepTuples)) *
 *              \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to do not give a
 *            non-empty range of increasing indices.
 *  \throw If \a a->getNbOfElems() != \a nbTargetValues.
 *  \throw If \a strictCompoCompare == \a true && \a a->getNumberOfComponents() !=
 *            \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *
 *  \ref cpp_mcdataarraydouble_setpartofvalues1 "Here is a Python example".
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
  const double *srcPt=a->getConstPointer();
  double *pt=getPointer()+bgTuples*nbComp+bgComp;
  if(assignTech)
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        for(int j=0;j<newNbOfComp;j++,srcPt++)
          pt[j*stepComp]=*srcPt;
    }
  else
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const double *srcPt2=srcPt;
          for(int j=0;j<newNbOfComp;j++,srcPt2++)
            pt[j*stepComp]=*srcPt2;
        }
    }
}

/*!
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step)..
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref cpp_mcdataarraydouble_setpartofvaluessimple1 "Here is a Python example".
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
 * Copy all values from another DataArrayDouble (\a a) into specified tuples and 
 * components of \a this array. Textual data is not copied.
 * The tuples and components to assign to are defined by C arrays of indices.
 * There are two *modes of usage*:
 * - If \a a->getNbOfElems() equals to number of values to assign to, then every value
 *   of \a a is assigned to its own location within \a this array. 
 * - If \a a includes one tuple, then all values of \a a are assigned to the specified
 *   components of every specified tuple of \a this array. In this mode it is required
 *   that \a a->getNumberOfComponents() equals to the number of specified components.
 *
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index <em>(pi)</em> varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \param [in] strictCompoCompare - this parameter is checked only if the
 *               *mode of usage* is the first; if it is \a true (default), 
 *               then \a a->getNumberOfComponents() must be equal 
 *               to the number of specified columns, else this is not required.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple/component given by <em>bgTuples / bgComp</em> is
 *         out of a valid range for \a this array.
 *  \throw In the first *mode of usage*, if <em>strictCompoCompare == true </em> and
 *         if <em> a->getNumberOfComponents() != (endComp - bgComp) </em>.
 *  \throw In the second *mode of usage*, if \a a->getNumberOfTuples() != 1 or
 *         <em> a->getNumberOfComponents() != (endComp - bgComp)</em>.
 *
 *  \ref cpp_mcdataarraydouble_setpartofvalues2 "Here is a Python example".
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
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tuples and components to assign to are defined by C arrays of indices.
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index (\a pi) varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index (\a pi) varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple/component given by <em>bgTuples / bgComp</em> is
 *         out of a valid range for \a this array.
 *
 *  \ref cpp_mcdataarraydouble_setpartofvaluessimple2 "Here is a Python example".
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
 * Copy all values from another DataArrayDouble (\a a) into specified tuples and 
 * components of \a this array. Textual data is not copied.
 * The tuples to assign to are defined by a C array of indices.
 * The components to assign to are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 * There are two *modes of usage*:
 * - If \a a->getNbOfElems() equals to number of values to assign to, then every value
 *   of \a a is assigned to its own location within \a this array. 
 * - If \a a includes one tuple, then all values of \a a are assigned to the specified
 *   components of every specified tuple of \a this array. In this mode it is required
 *   that \a a->getNumberOfComponents() equals to the number of specified components.
 *
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \param [in] strictCompoCompare - this parameter is checked only in the first
 *               *mode of usage*; if \a strictCompoCompare is \a true (default), 
 *               then \a a->getNumberOfComponents() must be equal 
 *               to the number of specified columns, else this is not required.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple given by \a bgTuples is out of a valid range for 
 *         \a this array.
 *  \throw In the first *mode of usage*, if <em>strictCompoCompare == true </em> and
 *         if <em> a->getNumberOfComponents()</em> is unequal to the number of components
 *         defined by <em>(bgComp,endComp,stepComp)</em>.
 *  \throw In the second *mode of usage*, if \a a->getNumberOfTuples() != 1 or
 *         <em> a->getNumberOfComponents()</em> is unequal to the number of components
 *         defined by <em>(bgComp,endComp,stepComp)</em>.
 *  \throw If parameters specifying components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref cpp_mcdataarraydouble_setpartofvalues3 "Here is a Python example".
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
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tuples to assign to are defined by a C array of indices.
 * The components to assign to are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple given by \a bgTuples is out of a valid range for 
 *         \a this array.
 *  \throw If parameters specifying components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref cpp_mcdataarraydouble_setpartofvaluessimple3 "Here is a Python example".
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
 * Copy all values from another DataArrayDouble into specified tuples and components
 * of \a this array. Textual data is not copied.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign values to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index (\a pi) varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \param [in] strictCompoCompare - if \a true (by default), then \a a->getNumberOfComponents() 
 *              must be equal to the number of columns to assign to, else an
 *              exception is thrown; if \a false, then it is only required that \a
 *              a->getNbOfElems() equals to number of values to assign to (this condition
 *              must be respected even if \a strictCompoCompare is \a true). The number of 
 *              values to assign to is given by following Python expression:
 *              \a nbTargetValues = 
 *              \c len(\c range(\a bgTuples,\a endTuples,\a stepTuples)) *
 *              \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to do not give a
 *            non-empty range of increasing indices.
 *  \throw If \a a->getNbOfElems() != \a nbTargetValues.
 *  \throw If \a strictCompoCompare == \a true && \a a->getNumberOfComponents() !=
 *            \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *
 */
void DataArrayDouble::setPartOfValues4(const DataArrayDouble *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValues4 : input DataArrayDouble is NULL !");
  const char msg[]="DataArrayDouble::setPartOfValues4";
  checkAllocated();
  a->checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=(int)std::distance(bgComp,endComp);
  int nbComp=getNumberOfComponents();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
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
  const double *srcPt=a->getConstPointer();
  double *pt=getPointer()+bgTuples*nbComp;
  if(assignTech)
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        for(const int *z=bgComp;z!=endComp;z++,srcPt++)
          pt[*z]=*srcPt;
    }
  else
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const double *srcPt2=srcPt;
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            pt[*z]=*srcPt2;
        }
    }
}

void DataArrayDouble::setPartOfValuesSimple4(double a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="DataArrayDouble::setPartOfValuesSimple4";
  checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int nbComp=getNumberOfComponents();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  double *pt=getPointer()+bgTuples*nbComp;
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(const int *z=bgComp;z!=endComp;z++)
      pt[*z]=a;
}

/*!
 * Copy some tuples from another DataArrayDouble into specified tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * Both the tuples to assign and the tuples to assign to are defined by a DataArrayInt.
 * All components of selected tuples are copied.
 *  \param [in] a - the array to copy values from.
 *  \param [in] tuplesSelec - the array specifying both source tuples of \a a and
 *              target tuples of \a this. \a tuplesSelec has two components, and the
 *              first component specifies index of the source tuple and the second
 *              one specifies index of the target tuple.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a tuplesSelec is NULL.
 *  \throw If \a tuplesSelec is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If \a tuplesSelec->getNumberOfComponents() != 2.
 *  \throw If any tuple index given by \a tuplesSelec is out of a valid range for 
 *         the corresponding (\a this or \a a) array.
 */
void DataArrayDouble::setPartOfValuesAdv(const DataArrayDouble *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception)
{
  if(!a || !tuplesSelec)
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
 * Copy some tuples from another DataArrayDouble (\a a) into contiguous tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * The tuples to assign to are defined by index of the first tuple, and
 * their number is defined by \a tuplesSelec->getNumberOfTuples().
 * The tuples to copy are defined by values of a DataArrayInt.
 * All components of selected tuples are copied.
 *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
 *              values to.
 *  \param [in] a - the array to copy values from.
 *  \param [in] tuplesSelec - the array specifying tuples of \a a to copy.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a tuplesSelec is NULL.
 *  \throw If \a tuplesSelec is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If \a tuplesSelec->getNumberOfComponents() != 1.
 *  \throw If <em>tupleIdStart + tuplesSelec->getNumberOfTuples() > this->getNumberOfTuples().</em>
 *  \throw If any tuple index given by \a tuplesSelec is out of a valid range for 
 *         \a a array.
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
 * Copy some tuples from another DataArrayDouble (\a a) into contiguous tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * The tuples to copy are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 * The tuples to assign to are defined by index of the first tuple, and
 * their number is defined by number of tuples to copy.
 * All components of selected tuples are copied.
 *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
 *              values to.
 *  \param [in] a - the array to copy values from.
 *  \param [in] bg - index of the first tuple to copy of the array \a a.
 *  \param [in] end2 - index of the tuple of \a a before which the tuples to copy
 *              are located.
 *  \param [in] step - index increment to get index of the next tuple to copy.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If <em>tupleIdStart + len(range(bg,end2,step)) > this->getNumberOfTuples().</em>
 *  \throw If parameters specifying tuples to copy, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for the array \a a.
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
 * Returns a value located at specified tuple and component.
 * This method is equivalent to DataArrayDouble::getIJ() except that validity of
 * parameters is checked. So this method is safe but expensive if used to go through
 * all values of \a this.
 *  \param [in] tupleId - index of tuple of interest.
 *  \param [in] compoId - index of component of interest.
 *  \return double - value located by \a tupleId and \a compoId.
 *  \throw If \a this is not allocated.
 *  \throw If condition <em>( 0 <= tupleId < this->getNumberOfTuples() )</em> is violated.
 *  \throw If condition <em>( 0 <= compoId < this->getNumberOfComponents() )</em> is violated.
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
 * Returns the last value of \a this. 
 *  \return double - the last value of \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this->getNumberOfTuples() < 1.
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

/*!
 * Sets a C array to be used as raw data of \a this. The previously set info
 *  of components is retained and re-sized. 
 * For more info see \ref MEDCouplingArraySteps1.
 *  \param [in] array - the C array to be used as raw data of \a this.
 *  \param [in] ownership - if \a true, \a array will be deallocated at destruction of \a this.
 *  \param [in] type - specifies how to deallocate \a array. If \a type == ParaMEDMEM::CPP_DEALLOC,
 *                     \c delete [] \c array; will be called. If \a type == ParaMEDMEM::C_DEALLOC,
 *                     \c free(\c array ) will be called.
 *  \param [in] nbOfTuple - new number of tuples in \a this.
 *  \param [in] nbOfCompo - new number of components in \a this.
 */
void DataArrayDouble::useArray(const double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
}

void DataArrayDouble::useExternalArrayWithRWAccess(const double *array, int nbOfTuple, int nbOfCompo)
{
  _info_on_compo.resize(nbOfCompo);
  _mem.useExternalArrayWithRWAccess(array,nbOfTuple*nbOfCompo);
  declareAsNew();
}

/*!
 * Checks if 0.0 value is present in \a this array. If it is the case, an exception
 * is thrown.
 * \throw If zero is found in \a this array.
 */
void DataArrayDouble::checkNoNullValues() const throw(INTERP_KERNEL::Exception)
{
  const double *tmp=getConstPointer();
  int nbOfElems=getNbOfElems();
  const double *where=std::find(tmp,tmp+nbOfElems,0.);
  if(where!=tmp+nbOfElems)
    throw INTERP_KERNEL::Exception("A value 0.0 have been detected !");
}

/*!
 * Computes minimal and maximal value in each component. An output array is filled
 * with \c 2 * \a this->getNumberOfComponents() values, so the caller is to allocate
 * enough memory before calling this method.
 *  \param [out] bounds - array of size at least 2 *\a this->getNumberOfComponents().
 *               It is filled as follows:<br>
 *               \a bounds[0] = \c min_of_component_0 <br>
 *               \a bounds[1] = \c max_of_component_0 <br>
 *               \a bounds[2] = \c min_of_component_1 <br>
 *               \a bounds[3] = \c max_of_component_1 <br>
 *               ...
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
  return bbox.retn();
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
void DataArrayDouble::computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayInt *& c, DataArrayInt *& cI) const throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeTupleIdsNearTuples : input pointer other is null !");
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> bbox=computeBBoxPerTuple(eps);
  other->checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  int otherNbOfCompo=other->getNumberOfComponents();
  if(nbOfCompo!=otherNbOfCompo)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeTupleIdsNearTuples : number of components should be equal between this and other !");
  int nbOfTuplesOther=other->getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cArr(DataArrayInt::New()),cIArr(DataArrayInt::New()); cArr->alloc(0,1); cIArr->pushBackSilent(0);
  switch(nbOfCompo)
    {
    case 3:
      {
        BBTree<3,int> myTree(bbox->getConstPointer(),0,0,getNumberOfTuples(),eps/10);
        FindTupleIdsNearTuplesAlg<3>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,cArr,cIArr);
        break;
      }
    case 2:
      {
        BBTree<2,int> myTree(bbox->getConstPointer(),0,0,getNumberOfTuples(),eps/10);
        FindTupleIdsNearTuplesAlg<2>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,cArr,cIArr);
        break;
      }
    case 1:
      {
        BBTree<1,int> myTree(bbox->getConstPointer(),0,0,getNumberOfTuples(),eps/10);
        FindTupleIdsNearTuplesAlg<1>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,cArr,cIArr);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Unexpected spacedim of coords for computeTupleIdsNearTuples. Must be 1, 2 or 3.");
    }
  c=cArr.retn(); cI=cIArr.retn();
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

/*!
 * Returns the maximal value and its location within \a this one-dimensional array.
 *  \param [out] tupleId - index of the tuple holding the maximal value.
 *  \return double - the maximal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Returns the maximal value within \a this array that is allowed to have more than
 *  one component.
 *  \return double - the maximal value among all values of \a this array.
 *  \throw If \a this is not allocated.
 */
double DataArrayDouble::getMaxValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *loc=std::max_element(begin(),end());
  return *loc;
}

/*!
 * Returns the maximal value and all its locations within \a this one-dimensional array.
 *  \param [out] tupleIds - a new instance of DataArrayInt containg indices of
 *               tuples holding the maximal value. The caller is to delete it using
 *               decrRef() as it is no more needed.
 *  \return double - the maximal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception)
{
  int tmp;
  tupleIds=0;
  double ret=getMaxValue(tmp);
  tupleIds=getIdsInRange(ret,ret);
  return ret;
}

/*!
 * Returns the minimal value and its location within \a this one-dimensional array.
 *  \param [out] tupleId - index of the tuple holding the minimal value.
 *  \return double - the minimal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Returns the minimal value within \a this array that is allowed to have more than
 *  one component.
 *  \return double - the minimal value among all values of \a this array.
 *  \throw If \a this is not allocated.
 */
double DataArrayDouble::getMinValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *loc=std::min_element(begin(),end());
  return *loc;
}

/*!
 * Returns the minimal value and all its locations within \a this one-dimensional array.
 *  \param [out] tupleIds - a new instance of DataArrayInt containg indices of
 *               tuples holding the minimal value. The caller is to delete it using
 *               decrRef() as it is no more needed.
 *  \return double - the minimal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception)
{
  int tmp;
  tupleIds=0;
  double ret=getMinValue(tmp);
  tupleIds=getIdsInRange(ret,ret);
  return ret;
}

/*!
 * Returns the average value of \a this one-dimensional array.
 *  \return double - the average value over all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
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

/*!
 * Returns the Euclidean norm of the vector defined by \a this array.
 *  \return double - the value of the Euclidean norm, i.e.
 *          the square root of the inner product of vector.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Returns the maximum norm of the vector defined by \a this array.
 *  \return double - the value of the maximum norm, i.e.
 *          the maximal absolute value among values of \a this array.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Accumulates values of each component of \a this array.
 *  \param [out] res - an array of length \a this->getNumberOfComponents(), allocated 
 *         by the caller, that is filled by this method with sum value for each
 *         component.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * This method returns the min distance from an external tuple defined by [ \a tupleBg , \a tupleEnd ) to \a this and
 * the first tuple in \a this that matches the returned distance. If there is no tuples in \a this an exception will be thrown.
 *
 *
 * \a this is expected to be allocated and expected to have a number of components equal to the distance from \a tupleBg to
 * \a tupleEnd. If not an exception will be thrown.
 *
 * \param [in] tupleBg start pointer (included) of input external tuple
 * \param [in] tupleEnd end pointer (not included) of input external tuple
 * \param [out] tupleId the tuple id in \a this that matches the min of distance between \a this and input external tuple
 * \return the min distance.
 * \sa MEDCouplingUMesh::distanceToPoint
 */
double DataArrayDouble::distanceToTuple(const double *tupleBg, const double *tupleEnd, int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbTuple=getNumberOfTuples();
  int nbComps=getNumberOfComponents();
  if(nbComps!=(int)std::distance(tupleBg,tupleEnd))
    { std::ostringstream oss; oss << "DataArrayDouble::distanceToTuple : size of input tuple is " << std::distance(tupleBg,tupleEnd) << " should be equal to the number of components in this : " << nbComps << " !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
  if(nbTuple==0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::distanceToTuple : no tuple in this ! No distance to compute !");
  double ret0=std::numeric_limits<double>::max();
  tupleId=-1;
  const double *work=getConstPointer();
  for(int i=0;i<nbTuple;i++)
    {
      double val=0.;
      for(int j=0;j<nbComps;j++,work++) 
        val+=(*work-tupleBg[j])*((*work-tupleBg[j]));
      if(val>=ret0)
        continue;
      else
        { ret0=val; tupleId=i; }
    }
  return sqrt(ret0);
}

/*!
 * Accumulate values of the given component of \a this array.
 *  \param [in] compId - the index of the component of interest.
 *  \return double - a sum value of \a compId-th component.
 *  \throw If \a this is not allocated.
 *  \throw If \a the condition ( 0 <= \a compId < \a this->getNumberOfComponents() ) is
 *         not respected.
 */
double DataArrayDouble::accumulate(int compId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *ptr=getConstPointer();
  int nbTuple=getNumberOfTuples();
  int nbComps=getNumberOfComponents();
  if(compId<0 || compId>=nbComps)
    throw INTERP_KERNEL::Exception("DataArrayDouble::accumulate : Invalid compId specified : No such nb of components !");
  double ret=0.;
  for(int i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

/*!
 * Converts each 2D point defined by the tuple of \a this array from the Polar to the
 * Cartesian coordinate system. The two components of the tuple of \a this array are 
 * considered to contain (1) radius and (2) angle of the point in the Polar CS.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          contains X and Y coordinates of the point in the Cartesian CS. The caller
 *          is to delete this array using decrRef() as it is no more needed. The array
 *          does not contain any textual info on components.
 *  \throw If \a this->getNumberOfComponents() != 2.
 */
DataArrayDouble *DataArrayDouble::fromPolarToCart() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Converts each 3D point defined by the tuple of \a this array from the Cylindrical to
 * the Cartesian coordinate system. The three components of the tuple of \a this array 
 * are considered to contain (1) radius, (2) azimuth and (3) altitude of the point in
 * the Cylindrical CS.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          contains X, Y and Z coordinates of the point in the Cartesian CS. The info
 *          on the third component is copied from \a this array. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this->getNumberOfComponents() != 3.
 */
DataArrayDouble *DataArrayDouble::fromCylToCart() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Converts each 3D point defined by the tuple of \a this array from the Spherical to
 * the Cartesian coordinate system. The three components of the tuple of \a this array 
 * are considered to contain (1) radius, (2) polar angle and (3) azimuthal angle of the
 * point in the Cylindrical CS.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          contains X, Y and Z coordinates of the point in the Cartesian CS. The info
 *          on the third component is copied from \a this array. The caller
 *          is to delete this array using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() != 3.
 */
DataArrayDouble *DataArrayDouble::fromSpherToCart() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes the doubly contracted product of every tensor defined by the tuple of \a this
 * array contating 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          is calculated from the tuple <em>(t)</em> of \a this array as follows:
 *          \f$ t[0]^2+t[1]^2+t[2]^2+2*t[3]^2+2*t[4]^2+2*t[5]^2\f$.
 *         The caller is to delete this result array using decrRef() as it is no more needed. 
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::doublyContractedProduct() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes the determinant of every square matrix defined by the tuple of \a this
 * array, which contains either 4, 6 or 9 components. The case of 6 components
 * corresponds to that of the upper triangular matrix.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          is the determinant of matrix of the corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed. 
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 */
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

/*!
 * Computes 3 eigenvalues of every upper triangular matrix defined by the tuple of
 * \a this array, which contains 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing 3
 *          components, whose each tuple contains the eigenvalues of the matrix of
 *          corresponding tuple of \a this array. 
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed. 
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::eigenValues() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes 3 eigenvectors of every upper triangular matrix defined by the tuple of
 * \a this array, which contains 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing 9
 *          components, whose each tuple contains 3 eigenvectors of the matrix of
 *          corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::eigenVectors() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes the inverse matrix of every matrix defined by the tuple of \a this
 * array, which contains either 4, 6 or 9 components. The case of 6 components
 * corresponds to that of the upper triangular matrix.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of components as \a this one, whose each tuple is the inverse
 *          matrix of the matrix of corresponding tuple of \a this array. 
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed. 
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 */
DataArrayDouble *DataArrayDouble::inverse() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes the trace of every matrix defined by the tuple of \a this
 * array, which contains either 4, 6 or 9 components. The case of 6 components
 * corresponds to that of the upper triangular matrix.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing 
 *          1 component, whose each tuple is the trace of
 *          the matrix of corresponding tuple of \a this array. 
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed. 
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 */
DataArrayDouble *DataArrayDouble::trace() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes the stress deviator tensor of every stress tensor defined by the tuple of
 * \a this array, which contains 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of components and tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::deviator() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Computes the magnitude of every vector defined by the tuple of
 * \a this array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and one component.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Computes the maximal value within every tuple of \a this array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and one component.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 */
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
  return ret.retn();
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
  return ret.retn();
}

/*!
 * Sorts value within every tuple of \a this array.
 *  \param [in] asc - if \a true, the values are sorted in ascending order, else,
 *              in descending order.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Converts every value of \a this array to its absolute value.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::abs() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  std::transform(ptr,ptr+nbOfElems,ptr,std::ptr_fun<double,double>(fabs));
  declareAsNew();
}

/*!
 * Apply a liner function to a given component of \a this array, so that
 * an array element <em>(x)</em> becomes \f$ a * x + b \f$.
 *  \param [in] a - the first coefficient of the function.
 *  \param [in] b - the second coefficient of the function.
 *  \param [in] compoId - the index of component to modify.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Apply a liner function to all elements of \a this array, so that
 * an element _x_ becomes \f$ a * x + b \f$.
 *  \param [in] a - the first coefficient of the function.
 *  \param [in] b - the second coefficient of the function.
 *  \throw If \a this is not allocated.
 */
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
 * Modify all elements of \a this array, so that
 * an element _x_ becomes \f$ numerator / x \f$.
 *  \param [in] numerator - the numerator used to modify array elements.
 *  \throw If \a this is not allocated.
 *  \throw If there is an element equal to 0.0 in \a this array.
 *  \warning If an exception is thrown because of presence of 0.0 element in \a this 
 *           array, all elements processed before detection of the zero element remain
 *           modified.
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
 * Returns a full copy of \a this array except that sign of all elements is reversed.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples and component as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
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

/*!
 * Returns a new DataArrayDouble created from \a this one by applying \a
 * FunctionToEvaluate to every tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] func - the \a FunctionToEvaluate declared as 
 *              \c bool (*\a func)(\c const \c double *\a pos, \c double *\a res), 
 *              where \a pos points to the first component of a tuple of \a this array
 *              and \a res points to the first component of a tuple of the result array.
 *              Note that length (number of components) of \a pos can differ from
 *              that of \a res.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a func returns \a false.
 */
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
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc1.
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and \a nbOfComp components.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If computing \a func fails.
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

/*!
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc0.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples and components as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If computing \a func fails.
 */
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
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc2.
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a func contains vars that are not in \a this->getInfoOnComponent().
 *  \throw If computing \a func fails.
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
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc3.
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] varsOrder - sequence of vars defining their order.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a func contains vars not in \a varsOrder.
 *  \throw If computing \a func fails.
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

/*!
 * Returns a new DataArrayInt contating indices of tuples of \a this one-dimensional
 * array whose values are within a given range. Textual data is not copied.
 *  \param [in] vmin - a lowest acceptable value.
 *  \param [in] vmax - a greatest acceptable value.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 1
 *
 *  \ref cpp_mcdataarraydouble_getidsinrange "Here is a C++ example".
 *
 *  \ref py_mcdataarraydouble_getidsinrange "Here is a Python example".
 */
DataArrayInt *DataArrayDouble::getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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

/*!
 * Returns a new DataArrayDouble by concatenating two given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If both \a a1 and \a a2 are NULL.
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents().
 */
DataArrayDouble *DataArrayDouble::Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> tmp(2);
  tmp[0]=a1; tmp[1]=a2;
  return Aggregate(tmp);
}

/*!
 * Returns a new DataArrayDouble by concatenating all given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If getNumberOfComponents() of arrays within \a arr.
 */
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

/*!
 * Returns a new DataArrayDouble by aggregating two given arrays, so that (1) the number
 * of components in the result array is a sum of the number of components of given arrays
 * and (2) the number of tuples in the result array is same as that of each of given
 * arrays. In other words the i-th tuple of result array includes all components of
 * i-th tuples of all given arrays.
 * Number of tuples in the given arrays must be  the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If both \a a1 and \a a2 are NULL.
 *  \throw If any given array is not allocated.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 */
DataArrayDouble *DataArrayDouble::Meld(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> arr(2);
  arr[0]=a1; arr[1]=a2;
  return Meld(arr);
}

/*!
 * Returns a new DataArrayDouble by aggregating all given arrays, so that (1) the number
 * of components in the result array is a sum of the number of components of given arrays
 * and (2) the number of tuples in the result array is same as that of each of given
 * arrays. In other words the i-th tuple of result array includes all components of
 * i-th tuples of all given arrays.
 * Number of tuples in the given arrays must be  the same.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If any given array is not allocated.
 *  \throw If getNumberOfTuples() of arrays within \a arr is different.
 */
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

/*!
 * Returns a new DataArrayDouble containing a dot product of two given arrays, so that
 * the i-th tuple of the result array is a sum of products of j-th components of i-th
 * tuples of given arrays (\f$ a_i = \sum_{j=1}^n a1_j * a2_j \f$).
 * Info on components and name is copied from the first of the given arrays.
 * Number of tuples and components in the given arrays must be the same.
 *  \param [in] a1 - a given array.
 *  \param [in] a2 - another given array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If any given array is not allocated.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents()
 */
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

/*!
 * Returns a new DataArrayDouble containing a cross product of two given arrays, so that
 * the i-th tuple of the result array contains 3 components of a vector which is a cross
 * product of two vectors defined by the i-th tuples of given arrays.
 * Info on components is copied from the first of the given arrays.
 * Number of tuples in the given arrays must be the same.
 * Number of components in the given arrays must be 3.
 *  \param [in] a1 - a given array.
 *  \param [in] a2 - another given array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != 3
 *  \throw If \a a2->getNumberOfComponents() != 3
 */
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

/*!
 * Returns a new DataArrayDouble containing maximal values of two given arrays.
 * Info on components is copied from the first of the given arrays.
 * Number of tuples and components in the given arrays must be the same.
 *  \param [in] a1 - an array to compare values with another one.
 *  \param [in] a2 - another array to compare values with the first one.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents()
 */
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

/*!
 * Returns a new DataArrayDouble containing minimal values of two given arrays.
 * Info on components is copied from the first of the given arrays.
 * Number of tuples and components in the given arrays must be the same.
 *  \param [in] a1 - an array to compare values with another one.
 *  \param [in] a2 - another array to compare values with the first one.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents()
 */
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

/*!
 * Returns a new DataArrayDouble that is a sum of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a sum of the corresponding values of \a a1 and \a a2,
 *   i.e.: _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - an array to sum up.
 *  \param [in] a2 - another array to sum up.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 */
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
  return ret.retn();
}

/*!
 * Adds values of another DataArrayDouble to values of \a this one. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   \a other array is added to the corresponding value of \a this array, i.e.:
 *   _a_ [ i, j ] += _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] += _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] += _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to add to \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 */
void DataArrayDouble::addEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::addEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::addEqual  !";
  checkAllocated();
  other->checkAllocated();
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

/*!
 * Returns a new DataArrayDouble that is a subtraction of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a subtraction of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - an array to subtract from.
 *  \param [in] a2 - an array to subtract.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 */
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
          return ret.retn();
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
          return ret.retn();
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
      return ret.retn();
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Substract !");//will always throw an exception
      return 0;
    }
}

/*!
 * Subtract values of another DataArrayDouble from values of \a this one. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   \a other array is subtracted from the corresponding value of \a this array, i.e.:
 *   _a_ [ i, j ] -= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] -= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] -= _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to subtract from \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 */
void DataArrayDouble::substractEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::substractEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::substractEqual  !";
  checkAllocated();
  other->checkAllocated();
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

/*!
 * Returns a new DataArrayDouble that is a product of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a product of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - a factor array.
 *  \param [in] a2 - another factor array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 */
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
  return ret.retn();
}

/*!
 * Multiply values of another DataArrayDouble to values of \a this one. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   \a other array is multiplied to the corresponding value of \a this array, i.e.:
 *   _a_ [ i, j ] *= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] *= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] *= _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to multiply to \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 */
void DataArrayDouble::multiplyEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::multiplyEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::multiplyEqual !";
  checkAllocated();
  other->checkAllocated();
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

/*!
 * Returns a new DataArrayDouble that is a division of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a division of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - a numerator array.
 *  \param [in] a2 - a denominator array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 *  \warning No check of division by zero is performed!
 */
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
          return ret.retn();
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
          return ret.retn();
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
      return ret.retn();
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Divide !");//will always throw an exception
      return 0;
    }
}

/*!
 * Divide values of \a this array by values of another DataArrayDouble. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *    \a this array is divided by the corresponding value of \a other one, i.e.:
 *   _a_ [ i, j ] /= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] /= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] /= _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to divide \a this one by.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 *  \warning No check of division by zero is performed!
 */
void DataArrayDouble::divideEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::divideEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::divideEqual !";
  checkAllocated();
  other->checkAllocated();
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
  std::ostringstream oss; oss.precision(17); oss << "(";
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
      ret->useExternalArrayWithRWAccess(_pt,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayDoubleTuple::buildDADouble : unable to build a requested DataArrayDouble instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Returns a new instance of DataArrayInt. The caller is to delete this array
 * using decrRef() as it is no more needed. 
 */
DataArrayInt *DataArrayInt::New()
{
  return new DataArrayInt;
}

/*!
 * Checks if raw data is allocated. Read more on the raw data
 * in \ref MEDCouplingArrayBasicsTuplesAndCompo "DataArrays infos" for more information.
 *  \return bool - \a true if the raw data is allocated, \a false else.
 */
bool DataArrayInt::isAllocated() const
{
  return getConstPointer()!=0;
}

/*!
 * Checks if raw data is allocated and throws an exception if it is not the case.
 *  \throw If the raw data is not allocated.
 */
void DataArrayInt::checkAllocated() const throw(INTERP_KERNEL::Exception)
{
  if(!isAllocated())
    throw INTERP_KERNEL::Exception("DataArrayInt::checkAllocated : Array is defined but not allocated ! Call alloc or setValues method first !");
}

std::size_t DataArrayInt::getHeapMemorySize() const
{
  std::size_t sz=(std::size_t)_mem.getNbOfElemAllocated();
  sz*=sizeof(int);
  return DataArray::getHeapMemorySize()+sz;
}

/*!
 * Sets information on all components. This method can change number of components
 * at certain conditions; if the conditions are not respected, an exception is thrown.
 * The number of components can be changed provided that \a this is not allocated.
 *
 * To know more on format of the component information see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - a vector of component infos.
 *  \throw If \a this->getNumberOfComponents() != \a info.size() && \a this->isAllocated()
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
 * Returns the only one value in \a this, if and only if number of elements
 * (nb of tuples * nb of components) is equal to 1, and that \a this is allocated.
 *  \return double - the sole value stored in \a this array.
 *  \throw If at least one of conditions stated above is not fulfilled.
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
 * Returns an integer value characterizing \a this array, which is useful for a quick
 * comparison of many instances of DataArrayInt.
 *  \return int - the hash value.
 *  \throw If \a this is not allocated.
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
 * Checks the number of tuples.
 *  \return bool - \a true if getNumberOfTuples() == 0, \a false else.
 *  \throw If \a this is not allocated.
 */
bool DataArrayInt::empty() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  return getNumberOfTuples()==0;
}

/*!
 * Returns a full copy of \a this. For more info on copying data arrays see
 * \ref MEDCouplingArrayBasicsCopyDeep.
 *  \return DataArrayInt * - a new instance of DataArrayInt.
 */
DataArrayInt *DataArrayInt::deepCpy() const
{
  return new DataArrayInt(*this);
}

/*!
 * Returns either a \a deep or \a shallow copy of this array. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeep and \ref MEDCouplingArrayBasicsCopyShallow.
 *  \param [in] dCpy - if \a true, a deep copy is returned, else, a shallow one.
 *  \return DataArrayInt * - either a new instance of DataArrayInt (if \a dCpy
 *          == \a true) or \a this instance (if \a dCpy == \a false).
 */
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

/*!
 * Copies all the data from another DataArrayInt. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeepAssign.
 *  \param [in] other - another instance of DataArrayInt to copy data from.
 *  \throw If the \a other is not allocated.
 */
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

/*!
 * This method reserve nbOfElems elements in memory ( nbOfElems*4 bytes ) \b without impacting the number of tuples in \a this.
 * If \a this has already been allocated, this method checks that \a this has only one component. If not an INTERP_KERNEL::Exception will be thrown.
 * If \a this has not already been allocated, number of components is set to one.
 * This method allows to reduce number of reallocations on invokation of DataArrayInt::pushBackSilent and DataArrayInt::pushBackValsSilent on \a this.
 * 
 * \sa DataArrayInt::pack, DataArrayInt::pushBackSilent, DataArrayInt::pushBackValsSilent
 */
void DataArrayInt::reserve(int nbOfElems) throw(INTERP_KERNEL::Exception)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    {
      _mem.reserve(nbOfElems);
    }
  else if(nbCompo==0)
    {
      _mem.reserve(nbOfElems);
      _info_on_compo.resize(1);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayInt::reserve : not available for DataArrayInt with number of components different than 1 !");
}

/*!
 * This method adds at the end of \a this the single value \a val. This method do \b not update its time label to avoid useless incrementation
 * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
 *
 * \param [in] val the value to be added in \a this
 * \throw If \a this has already been allocated with number of components different from one.
 * \sa DataArrayInt::pushBackValsSilent
 */
void DataArrayInt::pushBackSilent(int val) throw(INTERP_KERNEL::Exception)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    _mem.pushBack(val);
  else if(nbCompo==0)
    {
      _info_on_compo.resize(1);
      _mem.pushBack(val);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayInt::pushBackSilent : not available for DataArrayInt with number of components different than 1 !");
}

/*!
 * This method adds at the end of \a this a serie of values [\c valsBg,\c valsEnd). This method do \b not update its time label to avoid useless incrementation
 * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
 *
 *  \param [in] valsBg - an array of values to push at the end of \this.
 *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
 *              the last value of \a valsBg is \a valsEnd[ -1 ].
 * \throw If \a this has already been allocated with number of components different from one.
 * \sa DataArrayInt::pushBackSilent
 */
void DataArrayInt::pushBackValsSilent(const int *valsBg, const int *valsEnd) throw(INTERP_KERNEL::Exception)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    _mem.insertAtTheEnd(valsBg,valsEnd);
  else if(nbCompo==0)
    {
      _info_on_compo.resize(1);
      _mem.insertAtTheEnd(valsBg,valsEnd);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayInt::pushBackValsSilent : not available for DataArrayInt with number of components different than 1 !");
}

/*!
 * This method returns silently ( without updating time label in \a this ) the last value, if any and suppress it.
 * \throw If \a this is already empty.
 * \throw If \a this has number of components different from one.
 */
int DataArrayInt::popBackSilent() throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()==1)
    return _mem.popBack();
  else
    throw INTERP_KERNEL::Exception("DataArrayInt::popBackSilent : not available for DataArrayInt with number of components different than 1 !");
}

/*!
 * This method \b do \b not modify content of \a this. It only modify its memory footprint if the allocated memory is to high regarding real data to store.
 *
 * \sa DataArrayInt::getHeapMemorySize, DataArrayInt::reserve
 */
void DataArrayInt::pack() const throw(INTERP_KERNEL::Exception)
{
  _mem.pack();
}

/*!
 * Allocates the raw data in memory. If exactly as same memory as needed already
 * allocated, it is not re-allocated.
 *  \param [in] nbOfTuple - number of tuples of data to allocate.
 *  \param [in] nbOfCompo - number of components of data to allocate.
 *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
 */
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

/*!
 * Allocates the raw data in memory. If the memory was already allocated, then it is
 * freed and re-allocated. See an example of this method use
 * \ref MEDCouplingArraySteps1WC "here".
 *  \param [in] nbOfTuple - number of tuples of data to allocate.
 *  \param [in] nbOfCompo - number of components of data to allocate.
 *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
 */
void DataArrayInt::alloc(int nbOfTuple, int nbOfCompo) throw(INTERP_KERNEL::Exception)
{
  if(nbOfTuple<0 || nbOfCompo<0)
    throw INTERP_KERNEL::Exception("DataArrayInt::alloc : request for negative length of data !");
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*nbOfTuple);
  declareAsNew();
}

/*!
 * Assign zero to all values in \a this array. To know more on filling arrays see
 * \ref MEDCouplingArrayFill.
 * \throw If \a this is not allocated.
 */
void DataArrayInt::fillWithZero() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(0);
  declareAsNew();
}

/*!
 * Assign \a val to all values in \a this array. To know more on filling arrays see
 * \ref MEDCouplingArrayFill.
 *  \param [in] val - the value to fill with.
 *  \throw If \a this is not allocated.
 */
void DataArrayInt::fillWithValue(int val) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(val);
  declareAsNew();
}

/*!
 * Set all values in \a this array so that the i-th element equals to \a init + i
 * (i starts from zero). To know more on filling arrays see \ref MEDCouplingArrayFill.
 *  \param [in] init - value to assign to the first element of array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Returns a textual and human readable representation of \a this instance of
 * DataArrayInt. This text is shown when a DataArrayInt is printed in Python.
 *  \return std::string - text describing \a this DataArrayInt.
 */
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
  checkAllocated();
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

void DataArrayInt::reprCppStream(const char *varName, std::ostream& stream) const
{
  int nbTuples=getNumberOfTuples(),nbComp=getNumberOfComponents();
  const int *data=getConstPointer();
  stream << "DataArrayInt *" << varName << "=DataArrayInt::New();" << std::endl;
  if(nbTuples*nbComp>=1)
    {
      stream << "const int " << varName << "Data[" << nbTuples*nbComp << "]={";
      std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<int>(stream,","));
      stream << data[nbTuples*nbComp-1] << "};" << std::endl;
      stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
    }
  else
    stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
  stream << varName << "->setName(\"" << getName() << "\");" << std::endl;
}

/*!
 * Modifies \a this one-dimensional array so that each value \a v = \a indArrBg[ \a v ],
 * i.e. a current value is used as in index to get a new value from \a indArrBg.
 *  \param [in] indArrBg - pointer to the first element of array of new values to assign
 *         to \a this array.
 *  \param [in] indArrEnd - specifies the end of the array \a indArrBg, so that
 *              the last value of \a indArrBg is \a indArrEnd[ -1 ].
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If any value of \a this can't be used as a valid index for 
 *         [\a indArrBg, \a indArrEnd).
 */
void DataArrayInt::transformWithIndArr(const int *indArrBg, const int *indArrEnd) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Computes distribution of values of \a this one-dimensional array between given value
 * ranges (casts). This method is typically useful for entity number spliting by types,
 * for example. 
 *  \param [in] arrBg - the array of ascending values defining the value ranges. The i-th
 *         value of \a arrBg (\a arrBg[ i ]) gives the lowest value of the i-th range,
 *         and the greatest value of the i-th range equals to \a arrBg[ i+1 ] - 1. \a
 *         arrBg containing \a n values defines \a n-1 ranges. The last value of \a arrBg
 *         should be more than every value in \a this array.
 *  \param [in] arrEnd - specifies the end of the array \a arrBg, so that
 *              the last value of \a arrBg is \a arrEnd[ -1 ].
 *  \param [out] castArr - a new instance of DataArrayInt, of same size as \a this array
 *         (same number of tuples and components), the caller is to delete 
 *         using decrRef() as it is no more needed.
 *         This array contains indices of ranges for every value of \a this array. I.e.
 *         the i-th value of \a castArr gives the index of range the i-th value of \a this
 *         belongs to. Or, in other words, this parameter contains for each tuple in \a
 *         this in which cast it holds.
 *  \param [out] rankInsideCast - a new instance of DataArrayInt, of same size as \a this
 *         array, the caller is to delete using decrRef() as it is no more needed.
 *         This array contains ranks of values of \a this array within ranges
 *         they belongs to. I.e. the i-th value of \a rankInsideCast gives the rank of
 *         the i-th value of \a this array within the \a castArr[ i ]-th range, to which
 *         the i-th value of \a this belongs to. Or, in other words, this param contains 
 *         for each tuple its rank inside its cast. The rank is computed as difference
 *         between the value and the lowest value of range.
 *  \param [out] castsPresent - a new instance of DataArrayInt, containing indices of
 *         ranges (casts) to which at least one value of \a this array belongs.
 *         Or, in other words, this param contains the casts that \a this contains.
 *         The caller is to delete this array using decrRef() as it is no more needed.
 *
 * \b Example: If \a this contains [6,5,0,3,2,7,8,1,4] and \a arrBg contains [0,4,9] then
 *            the output of this method will be : 
 * - \a castArr       : [1,1,0,0,0,1,1,0,1]
 * - \a rankInsideCast: [2,1,0,3,2,3,4,1,0]
 * - \a castsPresent  : [0,1]
 *
 * I.e. values of \a this array belong to 2 ranges: #0 and #1. Value 6 belongs to the
 * range #1 and its rank within this range is 2; etc.
 *
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a arrEnd - arrBg < 2.
 *  \throw If any value of \a this is not less than \a arrEnd[-1].
 *  \warning The values contained in \a arrBg should be sorted ascendently. No
 *           check of this is be done. If not, the result is not warranted. 
 * 
 */
void DataArrayInt::splitByValueRange(const int *arrBg, const int *arrEnd,
                                     DataArrayInt *& castArr, DataArrayInt *& rankInsideCast, DataArrayInt *& castsPresent) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("Call splitByValueRange  method on DataArrayInt with only one component, you can call 'rearrange' method before !");
  int nbOfTuples=getNumberOfTuples();
  std::size_t nbOfCast=std::distance(arrBg,arrEnd);
  if(nbOfCast<2)
    throw INTERP_KERNEL::Exception("DataArrayInt::splitByValueRange : The input array giving the cast range values should be of size >=2 !");
  nbOfCast--;
  const int *work=getConstPointer();
  typedef std::reverse_iterator<const int *> rintstart;
  rintstart bg(arrEnd);//OK no problem because size of 'arr' is greater or equal 2
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
  castArr=ret1.retn();
  rankInsideCast=ret2.retn();
  castsPresent=ret3.retn();
}

/*!
 * Creates a one-dimensional DataArrayInt (\a res) whose contents are computed from 
 * values of \a this (\a a) and the given (\a indArr) arrays as follows:
 * \a res[ \a indArr[ \a a[ i ]]] = i. I.e. for each value in place i \a v = \a a[ i ],
 * new value in place \a indArr[ \a v ] is i.
 *  \param [in] indArrBg - the array holding indices within the result array to assign
 *         indices of values of \a this array pointing to values of \a indArrBg.
 *  \param [in] indArrEnd - specifies the end of the array \a indArrBg, so that
 *              the last value of \a indArrBg is \a indArrEnd[ -1 ].
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If any value of \a this array is not a valid index for \a indArrBg array.
 *  \throw If any value of \a indArrBg is not a valid index for \a this array.
 */
DataArrayInt *DataArrayInt::transformWithIndArrR(const int *indArrBg, const int *indArrEnd) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
      if(*pt>=0 && *pt<nbElemsIn)
        {
          int pos=indArrBg[*pt];
          if(pos>=0 && pos<nbOfTuples)
            tmp[pos]=i;
          else
            {
              std::ostringstream oss; oss << "DataArrayInt::transformWithIndArrR : error on tuple #" << i << " value of new pos is " << pos << " ( indArrBg[" << *pt << "]) ! Should be in [0," << nbOfTuples << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::transformWithIndArrR : error on tuple #" << i << " value is " << *pt << " and indirectionnal array as a size equal to " << nbElemsIn << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret.retn();
}

/*!
 * Creates a one-dimensional DataArrayInt of given length, whose contents are computed
 * from values of \a this array, which is supposed to contain a renumbering map in 
 * "Old to New" mode. The result array contains a renumbering map in "New to Old" mode.
 * To know how to use the renumbering maps see \ref MEDCouplingArrayRenumbering.
 *  \param [in] newNbOfElem - the number of tuples in the result array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 * 
 *  \ref cpp_mcdataarrayint_invertarrayo2n2n2o "Here is a C++ example".
 *
 *  \ref py_mcdataarrayint_invertarrayo2n2n2o "Here is a Python example".
 */
DataArrayInt *DataArrayInt::invertArrayO2N2N2O(int newNbOfElem) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(newNbOfElem,1);
  int nbOfOldNodes=getNumberOfTuples();
  const int *old2New=getConstPointer();
  int *pt=ret->getPointer();
  for(int i=0;i!=nbOfOldNodes;i++)
    if(old2New[i]!=-1)
      pt[old2New[i]]=i;
  return ret.retn();
}

/*!
 * This method is similar to DataArrayInt::invertArrayO2N2N2O except that 
 * Example : If \a this contains [0,1,2,0,3,4,5,4,6,4] this method will return [0,1,2,4,5,6,8] whereas DataArrayInt::invertArrayO2N2N2O returns [3,1,2,4,9,6,8]
 */
DataArrayInt *DataArrayInt::invertArrayO2N2N2OBis(int newNbOfElem) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(newNbOfElem,1);
  int nbOfOldNodes=getNumberOfTuples();
  const int *old2New=getConstPointer();
  int *pt=ret->getPointer();
  for(int i=nbOfOldNodes-1;i>=0;i--)
    if(old2New[i]!=-1)
      pt[old2New[i]]=i;
  return ret.retn();
}

/*!
 * Creates a one-dimensional DataArrayInt of given length, whose contents are computed
 * from values of \a this array, which is supposed to contain a renumbering map in 
 * "New to Old" mode. The result array contains a renumbering map in "Old to New" mode.
 * To know how to use the renumbering maps see \ref MEDCouplingArrayRenumbering.
 *  \param [in] newNbOfElem - the number of tuples in the result array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 * 
 *  \ref cpp_mcdataarrayint_invertarrayn2o2o2n "Here is a C++ example".
 *
 *  \ref py_mcdataarrayint_invertarrayn2o2o2n "Here is a Python example".
 */
DataArrayInt *DataArrayInt::invertArrayN2O2O2N(int oldNbOfElem) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(oldNbOfElem,1);
  const int *new2Old=getConstPointer();
  int *pt=ret->getPointer();
  std::fill(pt,pt+oldNbOfElem,-1);
  int nbOfNewElems=getNumberOfTuples();
  for(int i=0;i<nbOfNewElems;i++)
    pt[new2Old[i]]=i;
  return ret.retn();
}

bool DataArrayInt::isEqualIfNotWhy(const DataArrayInt& other, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,0,reason);
}

/*!
 * Checks if \a this and another DataArrayInt are fully equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayInt to compare with \a this one.
 *  \return bool - \a true if the two arrays are equal, \a false else.
 */
bool DataArrayInt::isEqual(const DataArrayInt& other) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,tmp);
}

/*!
 * Checks if values of \a this and another DataArrayInt are equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayInt to compare with \a this one.
 *  \return bool - \a true if the values of two arrays are equal, \a false else.
 */
bool DataArrayInt::isEqualWithoutConsideringStr(const DataArrayInt& other) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,0,tmp);
}

/*!
 * Checks if values of \a this and another DataArrayInt are equal. Comparison is
 * performed on sorted value sequences.
 * For more info see\ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayInt to compare with \a this one.
 *  \return bool - \a true if the sorted values of two arrays are equal, \a false else.
 */
bool DataArrayInt::isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a=deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> b=other.deepCpy();
  a->sort();
  b->sort();
  return a->isEqualWithoutConsideringStr(*b);
}

/*!
 * Sorts values of the array.
 *  \param [in] asc - \a true means ascending order, \a false, descending.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
void DataArrayInt::sort(bool asc) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::sort : only supported with 'this' array with ONE component !");
  _mem.sort(asc);
}

/*!
 * Reverse the array values.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
void DataArrayInt::reverse() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::reverse : only supported with 'this' array with ONE component !");
  _mem.reverse();
}

/*!
 * Checks that \a this array is consistently **increasing** or **decreasing** in value.
 * If not an exception is thrown.
 *  \param [in] increasing - if \a true, the array values should be increasing.
 *  \throw If sequence of values is not strictly monotonic in agreement with \a
 *         increasing arg.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
void DataArrayInt::checkMonotonic(bool increasing) const throw(INTERP_KERNEL::Exception)
{
  if(!isMonotonic(increasing))
    {
      if (increasing)
        throw INTERP_KERNEL::Exception("DataArrayInt::checkMonotonic : 'this' is not INCREASING monotonic !");
      else
        throw INTERP_KERNEL::Exception("DataArrayInt::checkMonotonic : 'this' is not DECREASING monotonic !");
    }
}

/*!
 * Checks that \a this array is consistently **increasing** or **decreasing** in value.
 *  \param [in] increasing - if \a true, array values should be increasing.
 *  \return bool - \a true if values change in accordance with \a increasing arg.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
bool DataArrayInt::isMonotonic(bool increasing) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::isMonotonic : only supported with 'this' array with ONE component !");
  int nbOfElements=getNumberOfTuples();
  const int *ptr=getConstPointer();
  if(nbOfElements==0)
    return true;
  int ref=ptr[0];
  if(increasing)
    {
      for(int i=1;i<nbOfElements;i++)
        {
          if(ptr[i]>=ref)
            ref=ptr[i];
          else
            return false;
        }
    }
  else
    {
      for(int i=1;i<nbOfElements;i++)
        {
          if(ptr[i]<=ref)
            ref=ptr[i];
          else
            return false;
        }
    }
  return true;
}

/*!
 * This method check that array consistently INCREASING or DECREASING in value.
 */
bool DataArrayInt::isStrictlyMonotonic(bool increasing) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::isStrictlyMonotonic : only supported with 'this' array with ONE component !");
  int nbOfElements=getNumberOfTuples();
  const int *ptr=getConstPointer();
  if(nbOfElements==0)
    return true;
  int ref=ptr[0];
  if(increasing)
    {
      for(int i=1;i<nbOfElements;i++)
        {
          if(ptr[i]>ref)
            ref=ptr[i];
          else
            return false;
        }
    }
  else
    {
      for(int i=1;i<nbOfElements;i++)
        {
          if(ptr[i]<ref)
            ref=ptr[i];
          else
            return false;
        }
    }
  return true;
}

/*!
 * This method check that array consistently INCREASING or DECREASING in value.
 */
void DataArrayInt::checkStrictlyMonotonic(bool increasing) const throw(INTERP_KERNEL::Exception)
{
  if(!isStrictlyMonotonic(increasing))
    {
      if (increasing)
        throw INTERP_KERNEL::Exception("DataArrayInt::checkStrictlyMonotonic : 'this' is not strictly INCREASING monotonic !");
      else
        throw INTERP_KERNEL::Exception("DataArrayInt::checkStrictlyMonotonic : 'this' is not strictly DECREASING monotonic !");
    }
}

/*!
 * Creates a new one-dimensional DataArrayInt of the same size as \a this and a given
 * one-dimensional arrays that must be of the same length. The result array describes
 * correspondence between \a this and \a other arrays, so that 
 * <em> other.getIJ(i,0) == this->getIJ(ret->getIJ(i),0)</em>. If such a permutation is
 * not possible because some element in \a other is not in \a this, an exception is thrown.
 *  \param [in] other - an array to compute permutation to.
 *  \return DataArrayInt * - a new instance of DataArrayInt, which is a permutation array
 * from \a this to \a other. The caller is to delete this array using decrRef() as it is
 * no more needed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a other->getNumberOfComponents() != 1.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples().
 *  \throw If \a other includes a value which is not in \a this array.
 * 
 *  \ref cpp_mcdataarrayint_buildpermutationarr "Here is a C++ example".
 *
 *  \ref py_mcdataarrayint_buildpermutationarr "Here is a Python example".
 */
DataArrayInt *DataArrayInt::buildPermutationArr(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1 || other.getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::buildPermutationArr : 'this' and 'other' have to have exactly ONE component !");
  int nbTuple=getNumberOfTuples();
  other.checkAllocated();
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
  return ret.retn();
}

/*!
 * Sets a C array to be used as raw data of \a this. The previously set info
 *  of components is retained and re-sized. 
 * For more info see \ref MEDCouplingArraySteps1.
 *  \param [in] array - the C array to be used as raw data of \a this.
 *  \param [in] ownership - if \a true, \a array will be deallocated at destruction of \a this.
 *  \param [in] type - specifies how to deallocate \a array. If \a type == ParaMEDMEM::CPP_DEALLOC,
 *                     \c delete [] \c array; will be called. If \a type == ParaMEDMEM::C_DEALLOC,
 *                     \c free(\c array ) will be called.
 *  \param [in] nbOfTuple - new number of tuples in \a this.
 *  \param [in] nbOfCompo - new number of components in \a this.
 */
void DataArrayInt::useArray(const int *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
}

void DataArrayInt::useExternalArrayWithRWAccess(const int *array, int nbOfTuple, int nbOfCompo)
{
  _info_on_compo.resize(nbOfCompo);
  _mem.useExternalArrayWithRWAccess(array,nbOfTuple*nbOfCompo);
  declareAsNew();
}

/*!
 * Returns a new DataArrayInt holding the same values as \a this array but differently
 * arranged in memory. If \a this array holds 2 components of 3 values:
 * \f$ x_0,x_1,x_2,y_0,y_1,y_2 \f$, then the result array holds these values arranged
 * as follows: \f$ x_0,y_0,x_1,y_1,x_2,y_2 \f$.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \warning Do not confuse this method with transpose()!
 */
DataArrayInt *DataArrayInt::fromNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayInt::fromNoInterlace : Not defined array !");
  int *tab=_mem.fromNoInterlace(getNumberOfComponents());
  DataArrayInt *ret=DataArrayInt::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

/*!
 * Returns a new DataArrayInt holding the same values as \a this array but differently
 * arranged in memory. If \a this array holds 2 components of 3 values:
 * \f$ x_0,y_0,x_1,y_1,x_2,y_2 \f$, then the result array holds these values arranged
 * as follows: \f$ x_0,x_1,x_2,y_0,y_1,y_2 \f$.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \warning Do not confuse this method with transpose()!
 */
DataArrayInt *DataArrayInt::toNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayInt::toNoInterlace : Not defined array !");
  int *tab=_mem.toNoInterlace(getNumberOfComponents());
  DataArrayInt *ret=DataArrayInt::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

/*!
 * Permutes values of \a this array as required by \a old2New array. The values are
 * permuted so that \c new[ \a old2New[ i ]] = \c old[ i ]. Number of tuples remains
 * the same as in \this one.
 * If a permutation reduction is needed, substr() or selectByTupleId() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *     giving a new position for i-th old value.
 */
void DataArrayInt::renumberInPlace(const int *old2New)
{
  checkAllocated();
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

/*!
 * Permutes values of \a this array as required by \a new2Old array. The values are
 * permuted so that \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of tuples remains
 * the same as in \this one.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
 *     giving a previous position of i-th new value.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
void DataArrayInt::renumberInPlaceR(const int *new2Old)
{
  checkAllocated();
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
 * Returns a copy of \a this array with values permuted as required by \a old2New array.
 * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ].
 * Number of tuples in the result array remains the same as in \this one.
 * If a permutation reduction is needed, renumberAndReduce() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *          giving a new position for i-th old value.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 */
DataArrayInt *DataArrayInt::renumber(const int *old2New) const
{
  checkAllocated();
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

/*!
 * Returns a copy of \a this array with values permuted as required by \a new2Old array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of
 * tuples in the result array remains the same as in \this one.
 * If a permutation reduction is needed, substr() or selectByTupleId() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
 *     giving a previous position of i-th new value.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayInt *DataArrayInt::renumberR(const int *new2Old) const
{
  checkAllocated();
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
 * Returns a shorten and permuted copy of \a this array. The new DataArrayInt is
 * of size \a newNbOfTuple and it's values are permuted as required by \a old2New array.
 * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ] for all
 * \a old2New[ i ] >= 0. In other words every i-th tuple in \a this array, for which 
 * \a old2New[ i ] is negative, is missing from the result array.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *     giving a new position for i-th old tuple and giving negative position for
 *     for i-th old tuple that should be omitted.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayInt *DataArrayInt::renumberAndReduce(const int *old2New, int newNbOfTuple) const
{
  checkAllocated();
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
 * Returns a shorten and permuted copy of \a this array. The new DataArrayInt is
 * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
 * \a new2OldBg array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
 * This method is equivalent to renumberAndReduce() except that convention in input is
 * \c new2old and \b not \c old2new.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
 *              tuple index in \a this array to fill the i-th tuple in the new array.
 *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
 *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
 *              \a new2OldBg <= \a pi < \a new2OldEnd.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayInt *DataArrayInt::selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  int nbComp=getNumberOfComponents();
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  int *pt=ret->getPointer();
  const int *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten and permuted copy of \a this array. The new DataArrayInt is
 * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
 * \a new2OldBg array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
 * This method is equivalent to renumberAndReduce() except that convention in input is
 * \c new2old and \b not \c old2new.
 * This method is equivalent to selectByTupleId() except that it prevents coping data
 * from behind the end of \a this array.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
 *              tuple index in \a this array to fill the i-th tuple in the new array.
 *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
 *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
 *              \a new2OldBg <= \a pi < \a new2OldEnd.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a new2OldEnd - \a new2OldBg > \a this->getNumberOfTuples().
 */
DataArrayInt *DataArrayInt::selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
  return ret.retn();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayInt contains every
 * (\a bg + \c i * \a step)-th tuple of \a this array located before the \a end2-th
 * tuple. Indices of the selected tuples are the same as ones returned by the Python
 * command \c range( \a bg, \a end2, \a step ).
 * This method is equivalent to selectByTupleIdSafe() except that the input array is
 * not constructed explicitly.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] bg - index of the first tuple to copy from \a this array.
 *  \param [in] end2 - index of the tuple before which the tuples to copy are located.
 *  \param [in] step - index increment to get index of the next tuple to copy.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If (\a end2 < \a bg) or (\a step <= 0).
 *  \sa DataArrayInt::substr.
 */
DataArrayInt *DataArrayInt::selectByTupleId2(int bg, int end2, int step) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  int nbComp=getNumberOfComponents();
  int newNbOfTuples=GetNumberOfItemGivenBES(bg,end2,step,"DataArrayInt::selectByTupleId2 : ");
  ret->alloc(newNbOfTuples,nbComp);
  int *pt=ret->getPointer();
  const int *srcPt=getConstPointer()+bg*nbComp;
  for(int i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
    std::copy(srcPt,srcPt+nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayInt contains ranges
 * of tuples specified by \a ranges parameter.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] ranges - std::vector of std::pair's each of which defines a range
 *              of tuples in [\c begin,\c end) format.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a end < \a begin.
 *  \throw If \a end > \a this->getNumberOfTuples().
 *  \throw If \a this is not allocated.
 */
DataArrayInt *DataArrayInt::selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuplesThis=getNumberOfTuples();
  if(ranges.empty())
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
      ret->alloc(0,nbOfComp);
      ret->copyStringInfoFrom(*this);
      return ret.retn();
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
  return ret.retn();
}

/*!
 * Returns a new DataArrayInt containing a renumbering map in "Old to New" mode.
 * This map, if applied to \a this array, would make it sorted. For example, if
 * \a this array contents are [9,10,0,6,4,11,3,7] then the contents of the result array
 * are [5,6,0,3,2,7,1,4]; if this result array (\a res) is used as an argument in call
 * \a this->renumber(\a res) then the returned array contains [0,3,4,6,7,9,10,11].
 * This method is useful for renumbering (in MED file for example). For more info
 * on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If there are equal values in \a this array.
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
 * Returns two arrays describing a surjective mapping from \a this set of values (\a A)
 * onto a set of values of size \a targetNb (\a B). The surjective function is 
 * \a B[ \a A[ i ]] = i. That is to say that for each \a id in [0,\a targetNb), where \a
 * targetNb < \a this->getNumberOfTuples(), there exists at least one tupleId (\a tid) so
 * that <em> this->getIJ( tid, 0 ) == id</em>. <br>
 * The first of out arrays returns indices of elements of \a this array, grouped by their
 * place in the set \a B. The second out array is the index of the first one; it shows how
 * many elements of \a A are mapped into each element of \a B. <br>
 * For more info on
 * mapping and its usage in renumbering see \ref MEDCouplingArrayRenumbering. <br>
 * \b Example:
 * - \a this: [0,3,2,3,2,2,1,2]
 * - \a targetNb: 4
 * - \a arr:  [0,  6,  2,4,5,7,  1,3]
 * - \a arrI: [0,1,2,6,8]
 *
 * This result means: <br>
 * the element of \a B 0 encounters within \a A once (\a arrI[ 0+1 ] - \a arrI[ 0 ]) and
 * its index within \a A is 0 ( \a arr[ 0:1 ] == \a arr[ \a arrI[ 0 ] : \a arrI[ 0+1 ]]);<br>
 * the element of \a B 2 encounters within \a A 4 times (\a arrI[ 2+1 ] - \a arrI[ 2 ]) and
 * its indices within \a A are [2,4,5,7] ( \a arr[ 2:6 ] == \a arr[ \a arrI[ 2 ] : 
 * \a arrI[ 2+1 ]]); <br> etc.
 *  \param [in] targetNb - the size of the set \a B. \a targetNb must be equal or more
 *         than the maximal value of \a A.
 *  \param [out] arr - a new instance of DataArrayInt returning indices of
 *         elements of \a this, grouped by their place in the set \a B. The caller is to delete
 *         this array using decrRef() as it is no more needed.
 *  \param [out] arrI - a new instance of DataArrayInt returning size of groups of equal
 *         elements of \a this. The caller is to delete this array using decrRef() as it
 *         is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If any value in \a this is more or equal to \a targetNb.
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
      if(tmp2>=0 && tmp2<targetNb)
        tmp[tmp2].push_back(i);
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::changeSurjectiveFormat : At pos " << i << " presence of element " << tmp2 << " ! should be in [0," << targetNb << ") !";
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
  arr=ret.retn();
  arrI=retI.retn();
}


/*!
 * Returns a new DataArrayInt containing a renumbering map in "Old to New" mode computed
 * from a zip representation of a surjective format (returned e.g. by
 * \ref ParaMEDMEM::DataArrayDouble::findCommonTuples() "DataArrayDouble::findCommonTuples()"
 * for example). The result array minimizes the permutation. <br>
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering. <br>
 * \b Example: <br>
 * - \a nbOfOldTuples: 10 
 * - \a arr          : [0,3, 5,7,9]
 * - \a arrIBg       : [0,2,5]
 * - \a newNbOfTuples: 7
 * - result array    : [0,1,2,0,3,4,5,4,6,4]
 *
 *  \param [in] nbOfOldTuples - number of tuples in the initial array \a arr.
 *  \param [in] arr - the array of tuple indices grouped by \a arrIBg array.
 *  \param [in] arrIBg - the array dividing all indices stored in \a arr into groups of
 *         (indices of) equal values. Its every element (except the last one) points to
 *         the first element of a group of equal values.
 *  \param [in] arrIEnd - specifies the end of \a arrIBg, so that the last element of \a
 *          arrIBg is \a arrIEnd[ -1 ].
 *  \param [out] newNbOfTuples - number of tuples after surjection application.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If any value of \a arr breaks condition ( 0 <= \a arr[ i ] < \a nbOfOldTuples ).
 */
DataArrayInt *DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(int nbOfOldTuples, const int *arr, const int *arrIBg, const int *arrIEnd, int &newNbOfTuples) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfOldTuples,1);
  int *pt=ret->getPointer();
  std::fill(pt,pt+nbOfOldTuples,-1);
  int nbOfGrps=((int)std::distance(arrIBg,arrIEnd))-1;
  const int *cIPtr=arrIBg;
  for(int i=0;i<nbOfGrps;i++)
    pt[arr[cIPtr[i]]]=-(i+2);
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
                {
                  if(arr[j]>=0 && arr[j]<nbOfOldTuples)
                    pt[arr[j]]=newNb;
                  else
                    {
                      std::ostringstream oss; oss << "DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2 : With element #" << j << " value is " << arr[j] << " should be in [0," << nbOfOldTuples << ") !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
              newNb++;
            }
        }
    }
  newNbOfTuples=newNb;
  return ret.retn();
}

/*!
 * Returns a new DataArrayInt containing a renumbering map in "New to Old" mode,
 * which if applied to \a this array would make it sorted ascendingly.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering. <br>
 * \b Example: <br>
 * - \a this: [2,0,1,1,0,1,2,0,1,1,0,0]
 * - result: [10,0,5,6,1,7,11,2,8,9,3,4]
 * - after applying result to \a this: [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2] 
 *
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
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
  return ret.retn();
}

/*!
 * Checks if contents of \a this array are equal to that of an array filled with
 * iota(). This method is particularly useful for DataArrayInt instances that represent
 * a renumbering array to check the real need in renumbering. 
 *  \return bool - \a true if \a this array contents == \a range( \a this->getNumberOfTuples())
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
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

/*!
 * Checks if all values in \a this array are equal to \a val.
 *  \param [in] val - value to check equality of array values to.
 *  \return bool - \a true if all values are \a val.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1
 */
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

/*!
 * Creates a new DataArrayDouble and assigns all (textual and numerical) data of \a this
 * array to the new one.
 *  \return DataArrayDouble * - the new instance of DataArrayInt.
 */
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
 * Returns a shorten copy of \a this array. The new DataArrayInt contains all
 * tuples starting from the \a tupleIdBg-th tuple and including all tuples located before
 * the \a tupleIdEnd-th one. This methods has a similar behavior as std::string::substr().
 * This method is a specialization of selectByTupleId2().
 *  \param [in] tupleIdBg - index of the first tuple to copy from \a this array.
 *  \param [in] tupleIdEnd - index of the tuple before which the tuples to copy are located.
 *          If \a tupleIdEnd == -1, all the tuples till the end of \a this array are copied.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a tupleIdBg < 0.
 *  \throw If \a tupleIdBg > \a this->getNumberOfTuples().
    \throw If \a tupleIdEnd != -1 && \a tupleIdEnd < \a this->getNumberOfTuples().
 *  \sa DataArrayInt::selectByTupleId2
 */
DataArrayInt *DataArrayInt::substr(int tupleIdBg, int tupleIdEnd) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Changes the number of components within \a this array so that its raw data **does
 * not** change, instead splitting this data into tuples changes.
 *  \param [in] newNbOfComp - number of components for \a this array to have.
 *  \throw If \a this is not allocated
 *  \throw If getNbOfElems() % \a newNbOfCompo != 0.
 *  \warning This method erases all (name and unit) component info set before!
 */
void DataArrayInt::rearrange(int newNbOfCompo) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfElems=getNbOfElems();
  if(nbOfElems%newNbOfCompo!=0)
    throw INTERP_KERNEL::Exception("DataArrayInt::rearrange : nbOfElems%newNbOfCompo!=0 !");
  _info_on_compo.clear();
  _info_on_compo.resize(newNbOfCompo);
  declareAsNew();
}

/*!
 * Changes the number of components within \a this array to be equal to its number
 * of tuples, and inversely its number of tuples to become equal to its number of 
 * components. So that its raw data **does not** change, instead splitting this
 * data into tuples changes.
 *  \throw If \a this is not allocated.
 *  \warning This method erases all (name and unit) component info set before!
 *  \warning Do not confuse this method with fromNoInterlace() and toNoInterlace()!
 *  \sa rearrange()
 */
void DataArrayInt::transpose() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  rearrange(nbOfTuples);
}

/*!
 * Returns a shorten or extended copy of \a this array. If \a newNbOfComp is less
 * than \a this->getNumberOfComponents() then the result array is shorten as each tuple
 * is truncated to have \a newNbOfComp components, keeping first components. If \a
 * newNbOfComp is more than \a this->getNumberOfComponents() then the result array is
 * expanded as each tuple is populated with \a dftValue to have \a newNbOfComp
 * components.  
 *  \param [in] newNbOfComp - number of components for the new array to have.
 *  \param [in] dftValue - value assigned to new values added to the new array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
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

/*!
 * Changes number of tuples in the array. If the new number of tuples is smaller
 * than the current number the array is truncated, otherwise the array is extended.
 *  \param [in] nbOfTuples - new number of tuples. 
 *  \throw If \a this is not allocated.
 */
void DataArrayInt::reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.reAlloc(getNumberOfComponents()*nbOfTuples);
  declareAsNew();
}


/*!
 * Returns a copy of \a this array composed of selected components.
 * The new DataArrayInt has the same number of tuples but includes components
 * specified by \a compoIds parameter. So that getNbOfElems() of the result array
 * can be either less, same or more than \a this->getNbOfElems().
 *  \param [in] compoIds - sequence of zero based indices of components to include
 *              into the new array.
 *  \return DataArrayInt * - the new instance of DataArrayInt that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If a component index (\a i) is not valid: 
 *         \a i < 0 || \a i >= \a this->getNumberOfComponents().
 *
 *  \ref cpp_mcdataarrayint_keepselectedcomponents "Here is a Python example".
 */
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
  return ret.retn();
}

/*!
 * Appends components of another array to components of \a this one, tuple by tuple.
 * So that the number of tuples of \a this array remains the same and the number of 
 * components increases.
 *  \param [in] other - the DataArrayInt to append to \a this one.
 *  \throw If \a this is not allocated.
 *  \throw If \a this and \a other arrays have different number of tuples.
 *
 *  \ref cpp_mcdataarrayint_meldwith "Here is a C++ example".
 *
 *  \ref py_mcdataarrayint_meldwith "Here is a Python example".
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

/*!
 * Copy all components in a specified order from another DataArrayInt.
 * The specified components become the first ones in \a this array.
 * Both numerical and textual data is copied. The number of tuples in \a this and
 * the other array can be different.
 *  \param [in] a - the array to copy data from.
 *  \param [in] compoIds - sequence of zero based indices of components, data of which is
 *              to be copied.
 *  \throw If \a a is NULL.
 *  \throw If \a compoIds.size() != \a a->getNumberOfComponents().
 *  \throw If \a compoIds[i] < 0 or \a compoIds[i] > \a this->getNumberOfComponents().
 *
 *  \ref cpp_mcdataarrayint_setselectedcomponents "Here is a Python example".
 */
void DataArrayInt::setSelectedComponents(const DataArrayInt *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setSelectedComponents : input DataArrayInt is NULL !");
  checkAllocated();
  a->checkAllocated();
  copyPartOfStringInfoFrom2(compoIds,*a);
  std::size_t partOfCompoSz=compoIds.size();
  int nbOfCompo=getNumberOfComponents();
  int nbOfTuples=std::min(getNumberOfTuples(),a->getNumberOfTuples());
  const int *ac=a->getConstPointer();
  int *nc=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(std::size_t j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}

/*!
 * Copy all values from another DataArrayInt into specified tuples and components
 * of \a this array. Textual data is not copied.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign values to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - index of the first component of \a this array to assign values to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \param [in] strictCompoCompare - if \a true (by default), then \a a->getNumberOfComponents() 
 *              must be equal to the number of columns to assign to, else an
 *              exception is thrown; if \a false, then it is only required that \a
 *              a->getNbOfElems() equals to number of values to assign to (this condition
 *              must be respected even if \a strictCompoCompare is \a true). The number of 
 *              values to assign to is given by following Python expression:
 *              \a nbTargetValues = 
 *              \c len(\c range(\a bgTuples,\a endTuples,\a stepTuples)) *
 *              \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to do not give a
 *            non-empty range of increasing indices.
 *  \throw If \a a->getNbOfElems() != \a nbTargetValues.
 *  \throw If \a strictCompoCompare == \a true && \a a->getNumberOfComponents() !=
 *            \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *
 *  \ref cpp_mcdataarrayint_setpartofvalues1 "Here is a Python example".
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
  int *pt=getPointer()+bgTuples*nbComp+bgComp;
  const int *srcPt=a->getConstPointer();
  if(assignTech)
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        for(int j=0;j<newNbOfComp;j++,srcPt++)
          pt[j*stepComp]=*srcPt;
    }
  else
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const int *srcPt2=srcPt;
          for(int j=0;j<newNbOfComp;j++,srcPt2++)
            pt[j*stepComp]=*srcPt2;
        }
    }
}

/*!
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step)..
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref cpp_mcdataarrayint_setpartofvaluessimple1 "Here is a Python example".
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
 * Copy all values from another DataArrayInt (\a a) into specified tuples and 
 * components of \a this array. Textual data is not copied.
 * The tuples and components to assign to are defined by C arrays of indices.
 * There are two *modes of usage*:
 * - If \a a->getNbOfElems() equals to number of values to assign to, then every value
 *   of \a a is assigned to its own location within \a this array. 
 * - If \a a includes one tuple, then all values of \a a are assigned to the specified
 *   components of every specified tuple of \a this array. In this mode it is required
 *   that \a a->getNumberOfComponents() equals to the number of specified components.
 * 
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index <em>(pi)</em> varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \param [in] strictCompoCompare - this parameter is checked only if the
 *               *mode of usage* is the first; if it is \a true (default), 
 *               then \a a->getNumberOfComponents() must be equal 
 *               to the number of specified columns, else this is not required.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple/component given by <em>bgTuples / bgComp</em> is
 *         out of a valid range for \a this array.
 *  \throw In the first *mode of usage*, if <em>strictCompoCompare == true </em> and
 *         if <em> a->getNumberOfComponents() != (endComp - bgComp) </em>.
 *  \throw In the second *mode of usage*, if \a a->getNumberOfTuples() != 1 or
 *         <em> a->getNumberOfComponents() != (endComp - bgComp)</em>.
 *
 *  \ref cpp_mcdataarrayint_setpartofvalues2 "Here is a Python example".
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
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tuples and components to assign to are defined by C arrays of indices.
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index (\a pi) varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index (\a pi) varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple/component given by <em>bgTuples / bgComp</em> is
 *         out of a valid range for \a this array.
 *
 *  \ref cpp_mcdataarrayint_setpartofvaluessimple2 "Here is a Python example".
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
 * Copy all values from another DataArrayInt (\a a) into specified tuples and 
 * components of \a this array. Textual data is not copied.
 * The tuples to assign to are defined by a C array of indices.
 * The components to assign to are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 * There are two *modes of usage*:
 * - If \a a->getNbOfElems() equals to number of values to assign to, then every value
 *   of \a a is assigned to its own location within \a this array. 
 * - If \a a includes one tuple, then all values of \a a are assigned to the specified
 *   components of every specified tuple of \a this array. In this mode it is required
 *   that \a a->getNumberOfComponents() equals to the number of specified components.
 *
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \param [in] strictCompoCompare - this parameter is checked only in the first
 *               *mode of usage*; if \a strictCompoCompare is \a true (default), 
 *               then \a a->getNumberOfComponents() must be equal 
 *               to the number of specified columns, else this is not required.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple given by \a bgTuples is out of a valid range for 
 *         \a this array.
 *  \throw In the first *mode of usage*, if <em>strictCompoCompare == true </em> and
 *         if <em> a->getNumberOfComponents()</em> is unequal to the number of components
 *         defined by <em>(bgComp,endComp,stepComp)</em>.
 *  \throw In the second *mode of usage*, if \a a->getNumberOfTuples() != 1 or
 *         <em> a->getNumberOfComponents()</em> is unequal to the number of components
 *         defined by <em>(bgComp,endComp,stepComp)</em>.
 *  \throw If parameters specifying components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref cpp_mcdataarrayint_setpartofvalues3 "Here is a Python example".
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
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tuples to assign to are defined by a C array of indices.
 * The components to assign to are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple given by \a bgTuples is out of a valid range for 
 *         \a this array.
 *  \throw If parameters specifying components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref cpp_mcdataarrayint_setpartofvaluessimple3 "Here is a Python example".
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

void DataArrayInt::setPartOfValues4(const DataArrayInt *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValues4 : input DataArrayInt is NULL !");
  const char msg[]="DataArrayInt::setPartOfValues4";
  checkAllocated();
  a->checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=(int)std::distance(bgComp,endComp);
  int nbComp=getNumberOfComponents();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
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
  const int *srcPt=a->getConstPointer();
  int *pt=getPointer()+bgTuples*nbComp;
  if(assignTech)
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        for(const int *z=bgComp;z!=endComp;z++,srcPt++)
          pt[*z]=*srcPt;
    }
  else
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const int *srcPt2=srcPt;
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            pt[*z]=*srcPt2;
        }
    }
}

void DataArrayInt::setPartOfValuesSimple4(int a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="DataArrayInt::setPartOfValuesSimple4";
  checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int nbComp=getNumberOfComponents();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  int *pt=getPointer()+bgTuples*nbComp;
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(const int *z=bgComp;z!=endComp;z++)
      pt[*z]=a;
}

/*!
 * Copy some tuples from another DataArrayInt into specified tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * Both the tuples to assign and the tuples to assign to are defined by a DataArrayInt.
 * All components of selected tuples are copied.
 *  \param [in] a - the array to copy values from.
 *  \param [in] tuplesSelec - the array specifying both source tuples of \a a and
 *              target tuples of \a this. \a tuplesSelec has two components, and the
 *              first component specifies index of the source tuple and the second
 *              one specifies index of the target tuple.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a tuplesSelec is NULL.
 *  \throw If \a tuplesSelec is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If \a tuplesSelec->getNumberOfComponents() != 2.
 *  \throw If any tuple index given by \a tuplesSelec is out of a valid range for 
 *         the corresponding (\a this or \a a) array.
 */
void DataArrayInt::setPartOfValuesAdv(const DataArrayInt *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception)
{
  if(!a || !tuplesSelec)
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
 * Copy some tuples from another DataArrayInt (\a a) into contiguous tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * The tuples to assign to are defined by index of the first tuple, and
 * their number is defined by \a tuplesSelec->getNumberOfTuples().
 * The tuples to copy are defined by values of a DataArrayInt.
 * All components of selected tuples are copied.
 *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
 *              values to.
 *  \param [in] a - the array to copy values from.
 *  \param [in] tuplesSelec - the array specifying tuples of \a a to copy.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a tuplesSelec is NULL.
 *  \throw If \a tuplesSelec is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If \a tuplesSelec->getNumberOfComponents() != 1.
 *  \throw If <em>tupleIdStart + tuplesSelec->getNumberOfTuples() > this->getNumberOfTuples().</em>
 *  \throw If any tuple index given by \a tuplesSelec is out of a valid range for 
 *         \a a array.
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
 * Copy some tuples from another DataArrayInt (\a a) into contiguous tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * The tuples to copy are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 * The tuples to assign to are defined by index of the first tuple, and
 * their number is defined by number of tuples to copy.
 * All components of selected tuples are copied.
 *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
 *              values to.
 *  \param [in] a - the array to copy values from.
 *  \param [in] bg - index of the first tuple to copy of the array \a a.
 *  \param [in] end2 - index of the tuple of \a a before which the tuples to copy
 *              are located.
 *  \param [in] step - index increment to get index of the next tuple to copy.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If <em>tupleIdStart + len(range(bg,end2,step)) > this->getNumberOfTuples().</em>
 *  \throw If parameters specifying tuples to copy, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for the array \a a.
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
 * Returns a value located at specified tuple and component.
 * This method is equivalent to DataArrayInt::getIJ() except that validity of
 * parameters is checked. So this method is safe but expensive if used to go through
 * all values of \a this.
 *  \param [in] tupleId - index of tuple of interest.
 *  \param [in] compoId - index of component of interest.
 *  \return double - value located by \a tupleId and \a compoId.
 *  \throw If \a this is not allocated.
 *  \throw If condition <em>( 0 <= tupleId < this->getNumberOfTuples() )</em> is violated.
 *  \throw If condition <em>( 0 <= compoId < this->getNumberOfComponents() )</em> is violated.
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
 * Returns the last value of \a this. 
 *  \return double - the last value of \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this->getNumberOfTuples() < 1.
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

/*!
 * Assign pointer to one array to a pointer to another appay. Reference counter of
 * \a arrayToSet is incremented / decremented.
 *  \param [in] newArray - the pointer to array to assign to \a arrayToSet.
 *  \param [in,out] arrayToSet - the pointer to array to assign to.
 */
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

/*!
 * Creates a new DataArrayInt containing IDs (indices) of tuples holding value equal to a
 * given one.
 *  \param [in] val - the value to find within \a this.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayInt::getIdsEqual(int val) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsEqual : the array must have only one component, you can call 'rearrange' method before !");
  const int *cptr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr==val)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * Creates a new DataArrayInt containing IDs (indices) of tuples holding value \b not
 * equal to a given one. 
 *  \param [in] val - the value to ignore within \a this.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayInt::getIdsNotEqual(int val) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsNotEqual : the array must have only one component, you can call 'rearrange' method before !");
  const int *cptr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr!=val)
      ret->pushBackSilent(i);
  return ret.retn();
}


/*!
 * Assigns \a newValue to all elements holding \a oldValue within \a this
 * one-dimensional array.
 *  \param [in] oldValue - the value to replace.
 *  \param [in] newValue - the value to assign.
 *  \return int - number of replacements performed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
int DataArrayInt::changeValue(int oldValue, int newValue) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::changeValue : the array must have only one component, you can call 'rearrange' method before !");
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

/*!
 * Creates a new DataArrayInt containing IDs (indices) of tuples holding value equal to
 * one of given values.
 *  \param [in] valsBg - an array of values to find within \a this array.
 *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
 *              the last value of \a valsBg is \a valsEnd[ -1 ].
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayInt::getIdsEqualList(const int *valsBg, const int *valsEnd) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsEqualList : the array must have only one component, you can call 'rearrange' method before !");
  std::set<int> vals2(valsBg,valsEnd);
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(vals2.find(*cptr)!=vals2.end())
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * Creates a new DataArrayInt containing IDs (indices) of tuples holding values \b not
 * equal to any of given values.
 *  \param [in] valsBg - an array of values to ignore within \a this array.
 *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
 *              the last value of \a valsBg is \a valsEnd[ -1 ].
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayInt::getIdsNotEqualList(const int *valsBg, const int *valsEnd) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsNotEqualList : the array must have only one component, you can call 'rearrange' method before !");
  std::set<int> vals2(valsBg,valsEnd);
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(vals2.find(*cptr)==vals2.end())
      ret->pushBackSilent(i);
  return ret.retn();
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
  checkAllocated();
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
  checkAllocated();
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
  checkAllocated();
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
  checkAllocated();
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
 * Returns \a true if a given value is present within \a this one-dimensional array.
 *  \param [in] value - the value to find within \a this array.
 *  \return bool - \a true in case if \a value is present within \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \sa locateValue()
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

/*!
 * Accumulates values of each component of \a this array.
 *  \param [out] res - an array of length \a this->getNumberOfComponents(), allocated 
 *         by the caller, that is filled by this method with sum value for each
 *         component.
 *  \throw If \a this is not allocated.
 */
void DataArrayInt::accumulate(int *res) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const int *ptr=getConstPointer();
  int nbTuple=getNumberOfTuples();
  int nbComps=getNumberOfComponents();
  std::fill(res,res+nbComps,0);
  for(int i=0;i<nbTuple;i++)
    std::transform(ptr+i*nbComps,ptr+(i+1)*nbComps,res,res,std::plus<int>());
}

int DataArrayInt::accumulate(int compId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const int *ptr=getConstPointer();
  int nbTuple=getNumberOfTuples();
  int nbComps=getNumberOfComponents();
  if(compId<0 || compId>=nbComps)
    throw INTERP_KERNEL::Exception("DataArrayInt::accumulate : Invalid compId specified : No such nb of components !");
  int ret=0;
  for(int i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

/*!
 * Returns a new DataArrayInt by concatenating two given arrays, so that (1) the number
 * of tuples in the result array is <em> a1->getNumberOfTuples() + a2->getNumberOfTuples() -
 * offsetA2</em> and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * First \a offsetA2 tuples of \a a2 are skipped and thus are missing from the result array.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \param [in] offsetA2 - number of tuples of \a a2 to skip.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents().
 */
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

/*!
 * Returns a new DataArrayInt by concatenating all given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If getNumberOfComponents() of arrays within \a arr.
 */
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

/*!
 * Returns the maximal value and its location within \a this one-dimensional array.
 *  \param [out] tupleId - index of the tuple holding the maximal value.
 *  \return double - the maximal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
int DataArrayInt::getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Returns the maximal value within \a this array that is allowed to have more than
 *  one component.
 *  \return int - the maximal value among all values of \a this array.
 *  \throw If \a this is not allocated.
 */
int DataArrayInt::getMaxValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const int *loc=std::max_element(begin(),end());
  return *loc;
}

/*!
 * Returns the minimal value and its location within \a this one-dimensional array.
 *  \param [out] tupleId - index of the tuple holding the minimal value.
 *  \return int - the minimal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
int DataArrayInt::getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
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
 * Returns the minimal value within \a this array that is allowed to have more than
 *  one component.
 *  \return int - the minimal value among all values of \a this array.
 *  \throw If \a this is not allocated.
 */
int DataArrayInt::getMinValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const int *loc=std::min_element(begin(),end());
  return *loc;
}

/*!
 * Converts every value of \a this array to its absolute value.
 *  \throw If \a this is not allocated.
 */
void DataArrayInt::abs() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  std::transform(ptr,ptr+nbOfElems,ptr,std::ptr_fun<int,int>(std::abs));
  declareAsNew();
}

/*!
 * Apply a liner function to a given component of \a this array, so that
 * an array element <em>(x)</em> becomes \f$ a * x + b \f$.
 *  \param [in] a - the first coefficient of the function.
 *  \param [in] b - the second coefficient of the function.
 *  \param [in] compoId - the index of component to modify.
 *  \throw If \a this is not allocated.
 */
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

/*!
 * Apply a liner function to all elements of \a this array, so that
 * an element _x_ becomes \f$ a * x + b \f$.
 *  \param [in] a - the first coefficient of the function.
 *  \param [in] b - the second coefficient of the function.
 *  \throw If \a this is not allocated.
 */
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
 * Returns a full copy of \a this array except that sign of all elements is reversed.
 *  \return DataArrayInt * - the new instance of DataArrayInt containing the
 *          same number of tuples and component as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
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
 * Modify all elements of \a this array, so that
 * an element _x_ becomes \f$ numerator / x \f$.
 *  \param [in] numerator - the numerator used to modify array elements.
 *  \throw If \a this is not allocated.
 *  \throw If there is an element equal to 0 in \a this array.
 *  \warning If an exception is thrown because of presence of 0 element in \a this 
 *           array, all elements processed before detection of the zero element remain
 *           modified.
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

/*!
 * Modify all elements of \a this array, so that
 * an element _x_ becomes \f$ x / val \f$.
 *  \param [in] val - the denominator used to modify array elements.
 *  \throw If \a this is not allocated.
 *  \throw If \a val == 0.
 */
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

/*!
 * Modify all elements of \a this array, so that
 * an element _x_ becomes  <em> x % val </em>.
 *  \param [in] val - the divisor used to modify array elements.
 *  \throw If \a this is not allocated.
 *  \throw If \a val <= 0.
 */
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
  checkAllocated();
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
 * Modify all elements of \a this array, so that
 * an element _x_ becomes <em> val % x </em>.
 *  \param [in] val - the divident used to modify array elements.
 *  \throw If \a this is not allocated.
 *  \throw If there is an element equal to or less than 0 in \a this array.
 *  \warning If an exception is thrown because of presence of an element <= 0 in \a this 
 *           array, all elements processed before detection of the zero element remain
 *           modified.
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

/*!
 * Returns a new DataArrayInt by aggregating two given arrays, so that (1) the number
 * of components in the result array is a sum of the number of components of given arrays
 * and (2) the number of tuples in the result array is same as that of each of given
 * arrays. In other words the i-th tuple of result array includes all components of
 * i-th tuples of all given arrays.
 * Number of tuples in the given arrays must be the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If both \a a1 and \a a2 are NULL.
 *  \throw If any given array is not allocated.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 */
DataArrayInt *DataArrayInt::Meld(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> arr(2);
  arr[0]=a1; arr[1]=a2;
  return Meld(arr);
}

/*!
 * Returns a new DataArrayInt by aggregating all given arrays, so that (1) the number
 * of components in the result array is a sum of the number of components of given arrays
 * and (2) the number of tuples in the result array is same as that of each of given
 * arrays. In other words the i-th tuple of result array includes all components of
 * i-th tuples of all given arrays.
 * Number of tuples in the given arrays must be  the same.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If any given array is not allocated.
 *  \throw If getNumberOfTuples() of arrays within \a arr is different.
 */
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
 * Returns a new DataArrayInt which is a minimal partition of elements of \a groups.
 * The i-th item of the result array is an ID of a set of elements belonging to a
 * unique set of groups, which the i-th element is a part of. This set of elements
 * belonging to a unique set of groups is called \a family, so the result array contains
 * IDs of families each element belongs to.
 *
 * \b Example: if we have two groups of elements: \a group1 [0,4] and \a group2 [ 0,1,2 ],
 * then there are 3 families:
 * - \a family1 (with ID 1) contains element [0] belonging to ( \a group1 + \a group2 ),
 * - \a family2 (with ID 2) contains elements [4] belonging to ( \a group1 ),
 * - \a family3 (with ID 3) contains element [1,2] belonging to ( \a group2 ), <br>
 * and the result array contains IDs of families [ 1,3,3,0,2 ]. <br> Note a family ID 0 which
 * stands for the element #3 which is in none of groups.
 *
 *  \param [in] groups - sequence of groups of element IDs.
 *  \param [in] newNb - total number of elements; it must be more than max ID of element
 *         in \a groups.
 *  \param [out] fidsOfGroups - IDs of families the elements of each group belong to.
 *  \return DataArrayInt * - a new instance of DataArrayInt containing IDs of families
 *         each element with ID from range [0, \a newNb ) belongs to. The caller is to
 *         delete this array using decrRef() as it is no more needed.
 *  \throw If any element ID in \a groups violates condition ( 0 <= ID < \a newNb ).
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
  return ret.retn();
}

/*!
 * Returns a new DataArrayInt which contains all elements of given one-dimensional
 * not negative arrays. The result array does not contain any duplicates and its values
 * are sorted in ascending order.
 *  \param [in] arr - sequence of DataArrayInt's to unite.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *         array using decrRef() as it is no more needed.
 *  \throw If any \a arr[i] is not allocated.
 *  \throw If \a arr[i]->getNumberOfComponents() != 1.
 *  \throw If any value of \a arr[i] is negative.
 */
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

/*!
 * Returns a new DataArrayInt which contains elements present in each of given one-dimensional
 * not negative arrays. The result array does not contain any duplicates and its values
 * are sorted in ascending order.
 *  \param [in] arr - sequence of DataArrayInt's to intersect.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *         array using decrRef() as it is no more needed.
 *  \throw If any \a arr[i] is not allocated.
 *  \throw If \a arr[i]->getNumberOfComponents() != 1.
 *  \throw If any value of \a arr[i] < 0.
 */
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
        throw INTERP_KERNEL::Exception("DataArrayInt::BuildIntersection : only single component allowed !");
      int tmp1;
      valm=std::min((*it)->getMinValue(tmp1),valm);
    }
  if(valm<0)
    throw INTERP_KERNEL::Exception("DataArrayInt::BuildIntersection : a negative value has been detected !");
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

/*!
 * Returns a new DataArrayInt which contains a complement of elements of \a this
 * one-dimensional array. I.e. the result array contains all elements from the range [0,
 * \a nbOfElement) not present in \a this array.
 *  \param [in] nbOfElement - maximal size of the result array.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *         array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If any element \a x of \a this array violates condition ( 0 <= \a x < \a
 *         nbOfElement ).
 */
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

/*!
 * Returns a new DataArrayInt containing elements of \a this one-dimensional missing
 * from an \a other one-dimensional array.
 *  \param [in] other - a DataArrayInt containing elements not to include in the result array.
 *  \return DataArrayInt * - a new instance of DataArrayInt with one component. The
 *         caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If \a other is NULL.
 *  \throw If \a other is not allocated.
 *  \throw If \a other->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \sa DataArrayInt::buildSubstractionOptimized()
 */
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

/*!
 * \a this is expected to have one component and to be sorted ascendingly (as for \a other).
 * \a other is expected to be a part of \a this. If not DataArrayInt::buildSubstraction should be called instead.
 * 
 * \param [in] other an array with one component and expected to be sorted ascendingly.
 * \ret list of ids in \a this but not in \a other.
 * \sa DataArrayInt::buildSubstraction
 */
DataArrayInt *DataArrayInt::buildSubstractionOptimized(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception)
{
  static const char *MSG="DataArrayInt::buildSubstractionOptimized : only single component allowed !";
  if(!other) throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstractionOptimized : NULL input array !");
  checkAllocated(); other->checkAllocated();
  if(getNumberOfComponents()!=1) throw INTERP_KERNEL::Exception(MSG);
  if(other->getNumberOfComponents()!=1) throw INTERP_KERNEL::Exception(MSG);
  const int *pt1Bg(begin()),*pt1End(end()),*pt2Bg(other->begin()),*pt2End(other->end()),*work1(pt1Bg),*work2(pt2Bg);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  for(;work1!=pt1End;work1++)
    {
      if(work2!=pt2End && *work1==*work2)
        work2++;
      else
        ret->pushBackSilent(*work1);
    }
  return ret.retn();
}


/*!
 * Returns a new DataArrayInt which contains all elements of \a this and a given
 * one-dimensional not negative arrays. The result array does not contain any duplicates
 * and its values are sorted in ascending order.
 *  \param [in] other - an array to unite with \a this one.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *         array using decrRef() as it is no more needed.
 *  \throw If \a this or \a other is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a other->getNumberOfComponents() != 1.
 *  \throw If any value of \a this or \a other is negative.
 */
DataArrayInt *DataArrayInt::buildUnion(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *>arrs(2);
  arrs[0]=this; arrs[1]=other;
  return BuildUnion(arrs);
}


/*!
 * Returns a new DataArrayInt which contains elements present in both \a this and a given
 * one-dimensional not negative arrays. The result array does not contain any duplicates
 * and its values are sorted in ascending order.
 *  \param [in] other - an array to intersect with \a this one.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *         array using decrRef() as it is no more needed.
 *  \throw If \a this or \a other is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a other->getNumberOfComponents() != 1.
 *  \throw If any value of \a this or \a other is negative.
 */
DataArrayInt *DataArrayInt::buildIntersection(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *>arrs(2);
  arrs[0]=this; arrs[1]=other;
  return BuildIntersection(arrs);
}

/*!
 * This method can be applied on allocated with one component DataArrayInt instance.
 * This method is typically relevant for sorted arrays. All consecutive duplicated items in \a this will appear only once in returned DataArrayInt instance.
 * Example : if \a this contains [1,2,2,3,3,3,3,4,5,5,7,7,7,19] the returned array will contain [1,2,3,4,5,7,19]
 * 
 * \return a newly allocated array that contain the result of the unique operation applied on \a this.
 * \throw if \a this is not allocated or if \a this has not exactly one component.
 */
DataArrayInt *DataArrayInt::buildUnique() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::buildUnique : only single component allowed !");
  int nbOfTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=deepCpy();
  int *data=tmp->getPointer();
  int *last=std::unique(data,data+nbOfTuples);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(std::distance(data,last),1);
  std::copy(data,last,ret->getPointer());
  return ret.retn();
}

/*!
 * Returns a new DataArrayInt which contains size of every of groups described by \a this
 * "index" array. Such "index" array is returned for example by 
 * \ref ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity
 * "MEDCouplingUMesh::buildDescendingConnectivity" and
 * \ref ParaMEDMEM::MEDCouplingUMesh::getNodalConnectivityIndex
 * "MEDCouplingUMesh::getNodalConnectivityIndex" etc.
 *  \return DataArrayInt * - a new instance of DataArrayInt, whose number of tuples
 *          equals to \a this->getNumberOfComponents() - 1, and number of components is 1.
 *          The caller is to delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this->getNumberOfTuples() < 2.
 *
 *  \b Example: <br> 
 *         - this contains [1,3,6,7,7,9,15]
 *         - result array contains [2,3,1,0,2,6],
 *          where 2 = 3 - 1, 3 = 6 - 3, 1 = 7 - 6 etc.
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
 * Modifies \a this one-dimensional array so that value of each element \a x
 * of \a this array (\a a) is computed as \f$ x_i = \sum_{j=0}^{i-1} a[ j ] \f$.
 * Or: for each i>0 new[i]=new[i-1]+old[i-1] for i==0 new[i]=0. Number of tuples
 * and components remains the same.<br>
 * This method is useful for allToAllV in MPI with contiguous policy. This method
 * differs from computeOffsets2() in that the number of tuples is \b not changed by
 * this one.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *
 *  \b Example: <br>
 *          - Before \a this contains [3,5,1,2,0,8]
 *          - After \a this contains  [0,3,8,9,11,11]<br>
 *          Note that the last element 19 = 11 + 8 is missing because size of \a this
 *          array is retained and thus there is no space to store the last element.
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
 * Modifies \a this one-dimensional array so that value of each element \a x
 * of \a this array (\a a) is computed as \f$ x_i = \sum_{j=0}^{i-1} a[ j ] \f$.
 * Or: for each i>0 new[i]=new[i-1]+old[i-1] for i==0 new[i]=0. Number
 * components remains the same and number of tuples is inceamented by one.<br>
 * This method is useful for allToAllV in MPI with contiguous policy. This method
 * differs from computeOffsets() in that the number of tuples is changed by this one.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *
 *  \b Example: <br>
 *          - Before \a this contains [3,5,1,2,0,8]
 *          - After \a this contains  [0,3,8,9,11,11,19]<br>
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
 * Returns a new DataArrayInt whose contents is computed from that of \a this and \a
 * offsets arrays as follows. \a offsets is a one-dimensional array considered as an
 * "index" array of a "iota" array, thus, whose each element gives an index of a group
 * beginning within the "iota" array. And \a this is a one-dimensional array
 * considered as a selector of groups described by \a offsets to include into the result array.
 *  \throw If \a offsets is NULL.
 *  \throw If \a offsets is not allocated.
 *  \throw If \a offsets->getNumberOfComponents() != 1.
 *  \throw If \a offsets is not monotonically increasing.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If any element of \a this is not a valid index for \a offsets array.
 *
 *  \b Example: <br>
 *          - \a this: [0,2,3]
 *          - \a offsets: [0,3,6,10,14,20]
 *          - result array: [0,1,2,6,7,8,9,10,11,12,13] == <br>
 *            \c range(0,3) + \c range(6,10) + \c range(10,14) ==<br>
 *            \c range( \a offsets[ \a this[0] ], offsets[ \a this[0]+1 ]) + 
 *            \c range( \a offsets[ \a this[1] ], offsets[ \a this[1]+1 ]) + 
 *            \c range( \a offsets[ \a this[2] ], offsets[ \a this[2]+1 ])
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
  return ret.retn();
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
  return ret.retn();
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
  return ret.retn();
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
  return ret.retn();
}

/*!
 * This method returns all different values found in \a this. This method throws if \a this has not been allocated.
 * But the number of components can be different from one.
 * \return a newly allocated array (that should be dealt by the caller) containing different values in \a this.
 */
DataArrayInt *DataArrayInt::getDifferentValues() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  std::set<int> ret;
  ret.insert(begin(),end());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2=DataArrayInt::New(); ret2->alloc((int)ret.size(),1);
  std::copy(ret.begin(),ret.end(),ret2->getPointer());
  return ret2.retn();
}

/*!
 * This method is a refinement of DataArrayInt::getDifferentValues because it returns not only different values in \a this but also, for each of
 * them it tells which tuple id have this id.
 * This method works only on arrays with one component (if it is not the case call DataArrayInt::rearrange(1) ).
 * This method returns two arrays having same size.
 * The instances of DataArrayInt in the returned vector have be specially allocated and computed by this method. Each of them should be dealt by the caller of this method.
 * Example : if this is equal to [1,0,1,2,0,2,2,-3,2] -> differentIds=[-3,0,1,2] and returned array will be equal to [[7],[1,4],[0,2],[3,5,6,8]]
 */
std::vector<DataArrayInt *> DataArrayInt::partitionByDifferentValues(std::vector<int>& differentIds) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::partitionByDifferentValues : this should have only one component !");
  int id=0;
  std::map<int,int> m,m2,m3;
  for(const int *w=begin();w!=end();w++)
    m[*w]++;
  differentIds.resize(m.size());
  std::vector<DataArrayInt *> ret(m.size());
  std::vector<int *> retPtr(m.size());
  for(std::map<int,int>::const_iterator it=m.begin();it!=m.end();it++,id++)
    {
      m2[(*it).first]=id;
      ret[id]=DataArrayInt::New();
      ret[id]->alloc((*it).second,1);
      retPtr[id]=ret[id]->getPointer();
      differentIds[id]=(*it).first;
    }
  id=0;
  for(const int *w=begin();w!=end();w++,id++)
    {
      retPtr[m2[*w]][m3[*w]++]=id;
    }
  return ret;
}

/*!
 * Returns a new DataArrayInt that is a sum of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a sum of the corresponding values of \a a1 and \a a2,
 *   i.e.: _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - an array to sum up.
 *  \param [in] a2 - another array to sum up.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 */
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
  return ret.retn();
}

/*!
 * Adds values of another DataArrayInt to values of \a this one. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   \a other array is added to the corresponding value of \a this array, i.e.:
 *   _a_ [ i, j ] += _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] += _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] += _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to add to \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 */
void DataArrayInt::addEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::addEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::addEqual  !";
  checkAllocated(); other->checkAllocated();
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

/*!
 * Returns a new DataArrayInt that is a subtraction of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a subtraction of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - an array to subtract from.
 *  \param [in] a2 - an array to subtract.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 */
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
          return ret.retn();
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
          return ret.retn();
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
      return ret.retn();
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Substract !");//will always throw an exception
      return 0;
    }
}

/*!
 * Subtract values of another DataArrayInt from values of \a this one. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   \a other array is subtracted from the corresponding value of \a this array, i.e.:
 *   _a_ [ i, j ] -= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] -= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] -= _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to subtract from \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 */
void DataArrayInt::substractEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::substractEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::substractEqual  !";
  checkAllocated(); other->checkAllocated();
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

/*!
 * Returns a new DataArrayInt that is a product of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a product of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - a factor array.
 *  \param [in] a2 - another factor array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 */
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
  return ret.retn();
}


/*!
 * Multiply values of another DataArrayInt to values of \a this one. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   \a other array is multiplied to the corresponding value of \a this array, i.e.:
 *   _a_ [ i, j ] *= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] *= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] *= _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to multiply to \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 */
void DataArrayInt::multiplyEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::multiplyEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::multiplyEqual !";
  checkAllocated(); other->checkAllocated();
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


/*!
 * Returns a new DataArrayInt that is a division of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a division of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - a numerator array.
 *  \param [in] a2 - a denominator array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 *  \warning No check of division by zero is performed!
 */
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
          return ret.retn();
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
          return ret.retn();
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
      return ret.retn();
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Divide !");//will always throw an exception
      return 0;
    }
}

/*!
 * Divide values of \a this array by values of another DataArrayInt. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *    \a this array is divided by the corresponding value of \a other one, i.e.:
 *   _a_ [ i, j ] /= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] /= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] /= _a2_ [ 0, j ].
 *
 *  \param [in] other - an array to divide \a this one by.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 *  \warning No check of division by zero is performed!
 */
void DataArrayInt::divideEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::divideEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::divideEqual !";
  checkAllocated(); other->checkAllocated();
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


/*!
 * Returns a new DataArrayInt that is a modulus of two given arrays. There are 3
 * valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *   the result array (_a_) is a division of the corresponding values of \a a1 and
 *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] % _a2_ [ i, j ].
 * 2.  The arrays have same number of tuples and one array, say _a2_, has one
 *   component. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] % _a2_ [ i, 0 ].
 * 3.  The arrays have same number of components and one array, say _a2_, has one
 *   tuple. Then
 *   _a_ [ i, j ] = _a1_ [ i, j ] % _a2_ [ 0, j ].
 *
 * Info on components is copied either from the first array (in the first case) or from
 * the array with maximal number of elements (getNbOfElems()).
 *  \param [in] a1 - a dividend array.
 *  \param [in] a2 - a divisor array.
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
 *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
 *         none of them has number of tuples or components equal to 1.
 *  \warning No check of division by zero is performed!
 */
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
          return ret.retn();
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
          return ret.retn();
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
      return ret.retn();
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Modulus !");//will always throw an exception
      return 0;
    }
}

/*!
 * Modify \a this array so that each value becomes a modulus of division of this value by
 * a value of another DataArrayInt. There are 3 valid cases.
 * 1.  The arrays have same number of tuples and components. Then each value of
 *    \a this array is divided by the corresponding value of \a other one, i.e.:
 *   _a_ [ i, j ] %= _other_ [ i, j ].
 * 2.  The arrays have same number of tuples and \a other array has one component. Then
 *   _a_ [ i, j ] %= _other_ [ i, 0 ].
 * 3.  The arrays have same number of components and \a other array has one tuple. Then
 *   _a_ [ i, j ] %= _a2_ [ 0, j ].
 *
 *  \param [in] other - a divisor array.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
 *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
 *         \a other has number of both tuples and components not equal to 1.
 *  \warning No check of division by zero is performed!
 */
void DataArrayInt::modulusEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::modulusEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::modulusEqual !";
  checkAllocated(); other->checkAllocated();
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

/*!
 * Returns a C array which is a renumbering map in "Old to New" mode for the input array.
 * This map, if applied to \a start array, would make it sorted. For example, if
 * \a start array contents are [9,10,0,6,4,11,3,7] then the contents of the result array is
 * [5,6,0,3,2,7,1,4].
 *  \param [in] start - pointer to the first element of the array for which the
 *         permutation map is computed.
 *  \param [in] end - pointer specifying the end of the array \a start, so that
 *         the last value of \a start is \a end[ -1 ].
 *  \return int * - the result permutation array that the caller is to delete as it is no
 *         more needed.
 *  \throw If there are equal values in the input array.
 */
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

/*!
 * Returns a new DataArrayInt containing an arithmetic progression
 * that is equal to the sequence returned by Python \c range(\a begin,\a  end,\a  step )
 * function.
 *  \param [in] begin - the start value of the result sequence.
 *  \param [in] end - limiting value, so that every value of the result array is less than
 *              \a end.
 *  \param [in] step - specifies the increment or decrement.
 *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a step == 0.
 *  \throw If \a end < \a begin && \a step > 0.
 *  \throw If \a end > \a begin && \a step < 0.
 */
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
  return ret.retn();
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
      ret->useExternalArrayWithRWAccess(_pt,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayIntTuple::buildDAInt : unable to build a requested DataArrayInt instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}
