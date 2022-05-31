// Copyright (C) 2007-2022  CEA/DEN, EDF R&D
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

#include "CommInterface.hxx"

namespace MEDCoupling
{
  MPI_Datatype ParaTraits<double>::MPIDataType = MPI_DOUBLE;

  MPI_Datatype ParaTraits<Int32>::MPIDataType = MPI_INT;

  MPI_Datatype ParaTraits<Int64>::MPIDataType = MPI_LONG;

  void CommInterface::gatherArrays(MPI_Comm comm, int root, const DataArrayIdType *array, std::vector< MCAuto<DataArrayIdType> >& arraysOut) const
  {
    this->gatherArraysT2<mcIdType>(comm,root,array,arraysOut);
  }

  /*!
   * Generalized AllGather collective communication.
   * This method send input \a array to all procs.
   */
  void CommInterface::allGatherArrays(MPI_Comm comm, const DataArrayIdType *array, std::vector< MCAuto<DataArrayIdType> >& arraysOut) const
  {
    this->allGatherArraysT2<mcIdType>(comm,array,arraysOut);
  }

  /*!
   * Generalized AllGather collective communication.
   * This method send input \a array to all procs.
   */
  int CommInterface::allGatherArrays(MPI_Comm comm, const DataArrayIdType *array, std::unique_ptr<mcIdType[]>& result, std::unique_ptr<mcIdType[]>& resultIndex) const
  {
    return this->allGatherArraysT<mcIdType>(comm,array,result,resultIndex);
  }

  void CommInterface::allToAllArrays(MPI_Comm comm, const std::vector< MCAuto<DataArrayDouble> >& arrays, std::vector< MCAuto<DataArrayDouble> >& arraysOut) const
  {
    this->allToAllArraysT<double>(comm,arrays,arraysOut);
  }

  void CommInterface::allToAllArrays(MPI_Comm comm, const std::vector< MCAuto<DataArrayDouble> >& arrays, MCAuto<DataArrayDouble>& arraysOut) const
  {
    std::unique_ptr<mcIdType[]> notUsed1;
    mcIdType notUsed2;
    this->allToAllArraysT2<double>(comm,arrays,arraysOut,notUsed1,notUsed2);
  }

  /*!
  * Generalized AllToAll collective communication.
  */
  void CommInterface::allToAllArrays(MPI_Comm comm, const std::vector< MCAuto<DataArrayIdType> >& arrays, std::vector< MCAuto<DataArrayIdType> >& arraysOut) const
  {
    this->allToAllArraysT<mcIdType>(comm,arrays,arraysOut);
  }

}
