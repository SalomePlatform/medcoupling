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

#ifndef __OVERLAPELEMENTLOCATOR_HXX__
#define __OVERLAPELEMENTLOCATOR_HXX__

#include "InterpolationOptions.hxx"
#include "MEDCouplingNatureOfField.hxx"

#include <mpi.h>
#include <vector>
#include <set>

namespace ParaMEDMEM
{
  class ParaFIELD;
  class ProcessorGroup;
  class ParaSUPPORT;
  class InterpolationMatrix;
  class MEDCouplingPointSet;
  class DataArrayInt;

  class OverlapElementLocator : public INTERP_KERNEL::InterpolationOptions
  {
  public:
    OverlapElementLocator(const ParaFIELD *sourceField, const ParaFIELD *targetField, const ProcessorGroup& group);
    virtual ~OverlapElementLocator();
    const MPI_Comm *getCommunicator() const;
  private:
    void _computeBoundingBoxes();
    bool _intersectsBoundingBox(int i, int j) const;
  private:
    const ParaFIELD *_local_source_field;
    const ParaFIELD *_local_target_field;
    int _local_space_dim;
    MEDCouplingPointSet *_local_source_mesh;
    MEDCouplingPointSet *_local_target_mesh;
    std::vector<MEDCouplingPointSet*> _distant_cell_meshes;
    std::vector<MEDCouplingPointSet*> _distant_face_meshes;
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
