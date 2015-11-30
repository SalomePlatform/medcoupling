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
  
  class OverlapElementLocator : public INTERP_KERNEL::InterpolationOptions
  {
  public:
    OverlapElementLocator(const ParaFIELD *sourceField, const ParaFIELD *targetField, const ProcessorGroup& group);
    virtual ~OverlapElementLocator();
    const MPI_Comm *getCommunicator() const;
    void exchangeMeshes(OverlapInterpolationMatrix& matrix);
    std::vector< std::pair<int,int> > getToDoList() const { return _to_do_list; }
    std::vector< std::vector< int > > getProcsInInteraction() const { return _proc_pairs; }
    std::string getSourceMethod() const;
    std::string getTargetMethod() const;
    const MEDCouplingPointSet *getSourceMesh(int procId) const;
    const DataArrayInt *getSourceIds(int procId) const;
    const MEDCouplingPointSet *getTargetMesh(int procId) const;
    const DataArrayInt *getTargetIds(int procId) const;
  private:
    void computeBoundingBoxes();
    bool intersectsBoundingBox(int i, int j) const;
    void sendLocalMeshTo(int procId, bool sourceOrTarget, OverlapInterpolationMatrix& matrix) const;
    void receiveRemoteMesh(int procId, bool sourceOrTarget);
    void sendMesh(int procId, const MEDCouplingPointSet *mesh, const DataArrayInt *idsToSend) const;
    void receiveMesh(int procId, MEDCouplingPointSet* &mesh, DataArrayInt *&ids) const;
  private:
    const ParaFIELD *_local_source_field;
    const ParaFIELD *_local_target_field;
    int _local_space_dim;
    MEDCouplingPointSet *_local_source_mesh;
    MEDCouplingPointSet *_local_target_mesh;
    std::vector<MEDCouplingPointSet*> _distant_cell_meshes;
    std::vector<MEDCouplingPointSet*> _distant_face_meshes;
    //! of size _group.size(). Contains for each source proc i, the ids of proc j the targets interact with. This vector is common for all procs in _group. 
    std::vector< std::vector< int > > _proc_pairs;
    //! list of interpolations couple to be done
    std::vector< std::pair<int,int> > _to_do_list;
    std::vector< std::pair<int,bool> > _procs_to_send;
    std::map<std::pair<int,bool>, MEDCouplingAutoRefCountObjectPtr< MEDCouplingPointSet > > _remote_meshes;
    std::map<std::pair<int,bool>, MEDCouplingAutoRefCountObjectPtr< DataArrayInt > > _remote_elems;
    double* _domain_bounding_boxes;
    const ProcessorGroup& _group;
    std::vector<int> _distant_proc_ids;
    const MPI_Comm *_comm;
    //Attributes only used by lazy side
    //std::vector<double> _values_added;
    //std::vector< std::vector<int> > _ids_per_working_proc;
    //std::vector< std::vector<int> > _ids_per_working_proc3;
    //std::vector< std::vector<double> > _values_per_working_proc;
  };

}

#endif
