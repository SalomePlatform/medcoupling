// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingMemArray.txx"

using namespace MEDCoupling;

template class MEDCoupling::MemArray<float>;
template class MEDCoupling::DataArrayTemplate<float>;
template class MEDCoupling::DataArrayTemplateClassic<float>;
template class MEDCoupling::DataArrayTemplateFP<float>;
template class MEDCoupling::DataArrayIterator<float>;

DataArrayFloat *DataArrayFloat::New()
{
  return new DataArrayFloat;
}

DataArrayFloat *DataArrayFloat::deepCopy() const
{
  return new DataArrayFloat(*this);
}

void DataArrayFloat::reprCppStream(const std::string& varName, std::ostream& stream) const
{
  int nbTuples(getNumberOfTuples()),nbComp(getNumberOfComponents());
  const float *data(begin());
  stream.precision(7);
  stream << "DataArrayFloat *" << varName << "=DataArrayFloat::New();" << std::endl;
  if(nbTuples*nbComp>=1)
    {
      stream << "const float " << varName << "Data[" << nbTuples*nbComp << "]={";
      std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<float>(stream,","));
      stream << data[nbTuples*nbComp-1] << "};" << std::endl;
      stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
    }
  else
    stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
  stream << varName << "->setName(\"" << getName() << "\");" << std::endl;
}

void DataArrayFloat::reprQuickOverview(std::ostream& stream) const
{
  static const std::size_t MAX_NB_OF_BYTE_IN_REPR=300;
  stream << "DataArrayFloat C++ instance at " << this << ". ";
  if(isAllocated())
    {
      int nbOfCompo=(int)_info_on_compo.size();
      if(nbOfCompo>=1)
        {
          int nbOfTuples=getNumberOfTuples();
          stream << "Number of tuples : " << nbOfTuples << ". Number of components : " << nbOfCompo << "." << std::endl;
          reprQuickOverviewData(stream,MAX_NB_OF_BYTE_IN_REPR);
        }
      else
        stream << "Number of components : 0.";
    }
  else
    stream << "*** No data allocated ****";
}

void DataArrayFloat::reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const
{
  const float *data(begin());
  int nbOfTuples(getNumberOfTuples());
  int nbOfCompo=(int)_info_on_compo.size();
  std::ostringstream oss2; oss2 << "[";
  oss2.precision(7);
  std::string oss2Str(oss2.str());
  bool isFinished=true;
  for(int i=0;i<nbOfTuples && isFinished;i++)
    {
      if(nbOfCompo>1)
        {
          oss2 << "(";
          for(int j=0;j<nbOfCompo;j++,data++)
            {
              oss2 << *data;
              if(j!=nbOfCompo-1) oss2 << ", ";
            }
          oss2 << ")";
        }
      else
        oss2 << *data++;
      if(i!=nbOfTuples-1) oss2 << ", ";
      std::string oss3Str(oss2.str());
      if(oss3Str.length()<maxNbOfByteInRepr)
        oss2Str=oss3Str;
      else
        isFinished=false;
    }
  stream << oss2Str;
  if(!isFinished)
    stream << "... ";
  stream << "]";
}

bool DataArrayFloat::isEqualIfNotWhy(const DataArrayFloat& other, float prec, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,prec,reason);
}

bool DataArrayFloat::isEqual(const DataArrayFloat& other, float prec) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,prec,tmp);
}

bool DataArrayFloat::isEqualWithoutConsideringStr(const DataArrayFloat& other, float prec) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,prec,tmp);
}

DataArrayFloatIterator *DataArrayFloat::iterator()
{
  return new DataArrayFloatIterator(this);
}

DataArrayFloatIterator::DataArrayFloatIterator(DataArrayFloat *da):DataArrayIterator<float>(da)
{
}

DataArrayFloatTuple::DataArrayFloatTuple(float *pt, int nbOfComp):DataArrayTuple<float>(pt,nbOfComp)
{
}

std::string DataArrayFloatTuple::repr() const
{
  std::ostringstream oss; oss.precision(7); oss << "(";
  for(int i=0;i<_nb_of_compo-1;i++)
    oss << _pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

float DataArrayFloatTuple::floatValue() const
{
  return this->zeValue();
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a MEDCoupling::DataArrayFloat::decrRef.
 * This method performs \b no copy of data. The content is only referenced using MEDCoupling::DataArrayFloat::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayFloat *DataArrayFloatTuple::buildDAFloat(int nbOfTuples, int nbOfCompo) const
{
  return this->buildDA(nbOfTuples,nbOfCompo);
}
