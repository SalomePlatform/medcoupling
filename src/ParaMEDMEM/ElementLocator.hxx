// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __ELEMENTLOCATOR_HXX__
#define __ELEMENTLOCATOR_HXX__

#include "InterpolationOptions.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MCType.hxx"

#include <mpi.h>
#include <vector>
#include <set>

namespace MEDCoupling
{
  class ParaFIELD;
  class ProcessorGroup;
  class InterpolationMatrix;
  class MEDCouplingPointSet;
  class DataArrayInt;

  /*! Internal class, not part of the public API. Used in InterpolationMatrix.
   *
   */
  class ElementLocator : public INTERP_KERNEL::InterpolationOptions
  {
  public:
    ElementLocator(const ParaFIELD& sourceField, const ProcessorGroup& distant_group, const ProcessorGroup& local_group);

    virtual ~ElementLocator();
    void exchangeMesh(int idistantrank,
                      MEDCouplingPointSet*& target_mesh,
                      int*& distant_ids);
    void exchangeMethod(const std::string& sourceMeth, int idistantrank, std::string& targetMeth);
    const std::vector<int>& getDistantProcIds() const { return _distant_proc_ids; }
    const MPI_Comm *getCommunicator() const;
    NatureOfField getLocalNature() const;
    //! This method is used to informed if there is -1D mesh on distant_group side or on local_group side.
    bool isM1DCorr() const { return _is_m1d_corr; }
    //Working side methods
    void recvPolicyFromLazySideW(std::vector<int>& policy);
    void sendSumToLazySideW(const std::vector< std::vector<int> >& distantLocEltIds, const std::vector< std::vector<double> >& partialSumRelToDistantIds);
    void recvSumFromLazySideW(std::vector< std::vector<double> >& globalSumRelToDistantIds);
    void sendCandidatesForAddElementsW(const std::vector<int>& distantGlobIds);
    void recvAddElementsFromLazyProcsW(std::vector<std::vector<int> >& elementsToAdd);
    //
    void sendLocalIdsToLazyProcsW(const std::vector< std::vector<int> >& distantLocEltIds);
    void recvGlobalIdsFromLazyProcsW(const std::vector< std::vector<int> >& distantLocEltIds, std::vector< std::vector<int> >& globalIds);
    void recvCandidatesGlobalIdsFromLazyProcsW(std::vector< std::vector<int> >& globalIds);
    void sendPartialSumToLazyProcsW(const std::vector<int>& distantGlobIds, const std::vector<double>& sum);
    //Lazy side methods
    int sendPolicyToWorkingSideL();
    void recvFromWorkingSideL();
    void sendToWorkingSideL();
    //
    void recvLocalIdsFromWorkingSideL();
    void sendGlobalIdsToWorkingSideL();
    void sendCandidatesGlobalIdsToWorkingSideL();
    //
    void recvSumFromWorkingSideL();
    void recvCandidatesForAddElementsL();
    void sendAddElementsToWorkingSideL();
  private:
    void _computeBoundingBoxes();
    bool _intersectsBoundingBox(int irank);
    void _exchangeMesh(MEDCouplingPointSet* local_mesh, MEDCouplingPointSet*& distant_mesh,
                       int iproc_distant, const DataArrayInt* distant_ids_send,
                       int*& distant_ids_recv);
  private:
    const ParaFIELD&  _local_para_field ;
    MEDCouplingPointSet* _local_cell_mesh;
    int _local_cell_mesh_space_dim;
    bool _is_m1d_corr;
    MEDCouplingPointSet* _local_face_mesh;
    std::vector<MEDCouplingPointSet*> _distant_cell_meshes;
    std::vector<MEDCouplingPointSet*> _distant_face_meshes;
    double* _domain_bounding_boxes;
    const ProcessorGroup& _distant_group;
    const ProcessorGroup& _local_group;
    ProcessorGroup* _union_group;
    std::vector<int> _distant_proc_ids;
    const MPI_Comm *_comm;
    //Attributes only used by lazy side
    std::vector<double> _values_added;
    std::vector< std::vector<int> > _ids_per_working_proc;
    std::vector< std::vector<int> > _ids_per_working_proc3;
    std::vector< std::vector<double> > _values_per_working_proc;
  public:
    static const int CUMULATIVE_POLICY=3;
    static const int NO_POST_TREATMENT_POLICY=7;
  };

}

#endif
