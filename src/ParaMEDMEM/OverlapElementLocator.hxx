// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __OVERLAPELEMENTLOCATOR_HXX__
#define __OVERLAPELEMENTLOCATOR_HXX__

#include "InterpolationOptions.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <mpi.h>
#include <vector>
#include <map>
#include <set>

namespace ParaMEDMEM
{
  class ParaFIELD;
  class ProcessorGroup;
  class OverlapInterpolationMatrix;
  
  typedef std::pair<int,int>   ProcCouple;     // a couple of proc IDs, typically used to define a exchange betw 2 procs

  class OverlapElementLocator : public INTERP_KERNEL::InterpolationOptions
  {
  public:
    OverlapElementLocator(const ParaFIELD *sourceField, const ParaFIELD *targetField, const ProcessorGroup& group,double epsAbs);
    virtual ~OverlapElementLocator();
    const MPI_Comm *getCommunicator() const;
    void exchangeMeshes(OverlapInterpolationMatrix& matrix);
    std::vector< std::pair<int,int> > getToDoList() const { return _to_do_list; }
    std::vector< int > getProcsToSendFieldData() const { return _procs_to_send_field; }
    std::string getSourceMethod() const;
    std::string getTargetMethod() const;
    const MEDCouplingPointSet *getSourceMesh(int procId) const;
    const DataArrayInt *getSourceIds(int procId) const;
    const MEDCouplingPointSet *getTargetMesh(int procId) const;
    const DataArrayInt *getTargetIds(int procId) const;
  private:
    void computeBoundingBoxesAndTodoList();
    bool intersectsBoundingBox(int i, int j) const;
    void sendLocalMeshTo(int procId, bool sourceOrTarget, OverlapInterpolationMatrix& matrix) const;
    void receiveRemoteMeshFrom(int procId, bool sourceOrTarget);
    void sendMesh(int procId, const MEDCouplingPointSet *mesh, const DataArrayInt *idsToSend) const;
    void receiveMesh(int procId, MEDCouplingPointSet* &mesh, DataArrayInt *&ids) const;
  private:
    typedef MEDCouplingAutoRefCountObjectPtr< MEDCouplingPointSet >  AutoMCPointSet;
    typedef MEDCouplingAutoRefCountObjectPtr< DataArrayInt >         AutoDAInt;
    typedef std::pair<int,bool>  Proc_SrcOrTgt;  // a key indicating a proc ID and whether the data is for source mesh/field or target mesh/field

    static const int START_TAG_MESH_XCH;

    const ParaFIELD *_local_source_field;
    const ParaFIELD *_local_target_field;

    int _local_space_dim;
    MEDCouplingPointSet *_local_source_mesh;
    MEDCouplingPointSet *_local_target_mesh;

    /*! of size _group.size(). Contains for each source proc i, the ids of proc j the targets interact with.
        This vector is common for all procs in _group. */
    std::vector< std::vector< int > > _proc_pairs;
    //! list of interpolation couples to be done by this proc only. This is a simple extraction of the member _pairsToBeDonePerProc
    std::vector< ProcCouple > _to_do_list;
    //! list of procs the local proc will have to send mesh data to:
    std::vector< Proc_SrcOrTgt > _procs_to_send_mesh;
    /*! list of procs the local proc will have to send field data to for the final matrix-vector computation:
     * This can be different from _procs_to_send_mesh because interpolation matrix bits are computed on a potentially
     * different proc than the target one.   */
    std::vector< int > _procs_to_send_field;
    //! Set of distant meshes
    std::map< Proc_SrcOrTgt,  AutoMCPointSet > _remote_meshes;
    //! Set of cell ID mappings for the above distant meshes (because only part of the meshes are exchanged)
    std::map< Proc_SrcOrTgt, AutoDAInt > _remote_elems;
    double* _domain_bounding_boxes;
    //! bounding box absolute adjustment
    double _epsAbs;

    std::vector<int> _distant_proc_ids;

    const ProcessorGroup& _group;
    const MPI_Comm *_comm;
  };

}

#endif
