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
#ifndef __GLOBALIZERMESH_HXX__
#define __GLOBALIZERMESH_HXX__

#include "MEDCouplingNatureOfField.hxx"

#include <mpi.h>
#include <string>
#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingFieldDouble;

  class GlobalizerMesh
  {
  protected:
    GlobalizerMesh(const MPI_Comm *comm, MEDCouplingFieldDouble *localField);
  public:
    NatureOfField getLocalNature() const;
    virtual ~GlobalizerMesh();
  protected:
    const MPI_Comm *_comm;
    MEDCouplingFieldDouble *_local_field;
  };

  class GlobalizerMeshWorkingSide : public GlobalizerMesh
  {
  public:
    GlobalizerMeshWorkingSide(const MPI_Comm *comm, MEDCouplingFieldDouble *localField, const std::string& distantMeth, const std::vector<int>& lazyProcs);
    ~GlobalizerMeshWorkingSide();
    //
    const std::vector<int>& getProcIdsInInteraction() const;
    void sendSumToLazySide(const std::vector< std::vector<int> >& distantLocEltIds, const std::vector< std::vector<double> >& partialSumRelToDistantIds);
    void recvSumFromLazySide(std::vector< std::vector<double> >& globalSumRelToDistantIds);
  private:
    std::string _distant_method;
    std::vector<int> _lazy_procs;
  };

  class GlobalizerMeshLazySide : public GlobalizerMesh
  {
  public:
    GlobalizerMeshLazySide(const MPI_Comm *comm, MEDCouplingFieldDouble *localField, const std::vector<int>& computeProcs);
    ~GlobalizerMeshLazySide();
    void recvFromWorkingSide();
    void sendToWorkingSide();
  private:
    std::vector<int> _compute_procs;
    std::vector<double> _values_added;
    std::vector< std::vector<int> > _ids_per_working_proc;
  };
}

#endif
