// Copyright (C) 2020-2024  CEA, EDF
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

#pragma once

#include "ParaDataArray.hxx"
#include "CommInterface.hxx"
#include "MEDCouplingMemArray.txx"

#include <sstream>

namespace MEDCoupling
{
  template<class T>
  ParaDataArrayTemplate<T>::ParaDataArrayTemplate(typename Traits<T>::ArrayType *seqDa)
  {
    this->_seq_da.takeRef(seqDa);
  }

  template<class T>
  std::size_t ParaDataArrayTemplate<T>::getHeapMemorySizeWithoutChildren() const
  {
    return 0;
  }

  template<class T>
  std::vector<const BigMemoryObject *> ParaDataArrayTemplate<T>::getDirectChildrenWithNull() const
  {
    return { this->_seq_da };
  }
  
  template<class T>
  void ParaDataArrayTemplate<T>::checkOKOneComponent(const std::string& msg) const
  {
    if(this->_seq_da.isNull())
    {
      std::ostringstream oss; oss << msg << " : nullptr internal pointer !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
    this->_seq_da->checkAllocated();
    if( this->_seq_da->getNumberOfComponents()!=1 )
    {
      std::ostringstream oss; oss << msg << " : internal seq dataarray does not contain one component as expected !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  }

  /*!
    Parallel version of DataArrayInt::buildComplement. Returns result on proc 0. Not allocated DataArrayT is returned for all procs.
   */
  template<class T>
  DataArrayIdType *ParaDataArrayDiscrete<T>::buildComplement(T nbOfElems) const
  {
    using DataArrayT = typename Traits<T>::ArrayType;
    this->checkOKOneComponent("ParaDataArray::buildComplement");
    MPI_Comm comm(MPI_COMM_WORLD);
    CommInterface ci;
    int size;
    ci.commSize(comm,&size);
    std::vector< MCAuto<DataArrayT> > idsCaptured(size);
    for(int curRk = 0 ; curRk < size ; ++curRk)
    {
      T curStart(0),curEnd(0);
      DataArrayTools<T>::GetSlice(0,nbOfElems,1,ToIdType(curRk),ToIdType(size),curStart,curEnd);
      MCAuto<DataArrayIdType> idsInGlobalIds(this->_seq_da->findIdsInRange(curStart,curEnd));
      idsCaptured[curRk] = this->_seq_da->selectByTupleIdSafe(idsInGlobalIds->begin(),idsInGlobalIds->end());
    }
    // communication : 1 arrays are going to be all2allized : ids
    MCAuto<DataArrayT> aggregatedIds;
    {
      std::vector< MCAuto<DataArrayT> > myRkIdsCaptured;
      ci.allToAllArraysT<T>(comm,idsCaptured,myRkIdsCaptured);
      aggregatedIds = DataArrayT::Aggregate(FromVecAutoToVecOfConst<DataArrayT>(myRkIdsCaptured));
    }
    aggregatedIds->sort();
    aggregatedIds = aggregatedIds->buildUnique();
    int rank(-1);
    ci.commRank(comm,&rank);
    T vmin(std::numeric_limits<T>::max()),vmax(-std::numeric_limits<T>::max());
    DataArrayTools<T>::GetSlice(0,nbOfElems,1,ToIdType(rank),ToIdType(size),vmin,vmax);
    aggregatedIds->applyLin(1,-vmin);
    MCAuto<DataArrayIdType> seqComp(aggregatedIds->buildComplement(ToIdType(vmax-vmin)));
    seqComp->applyLin(1,ToIdType(vmin));
    //
    std::vector< MCAuto<DataArrayIdType> > arraysOut;
    ci.gatherArrays(comm,0,seqComp,arraysOut);
    MCAuto<DataArrayIdType> ret(DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(arraysOut)));
    return ret.retn();
  }

}
