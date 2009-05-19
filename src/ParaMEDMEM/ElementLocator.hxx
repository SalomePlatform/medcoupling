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
#ifndef __ELEMENTLOCATOR_HXX__
#define __ELEMENTLOCATOR_HXX__

#include "InterpolationOptions.hxx"

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

  class ElementLocator : public INTERP_KERNEL::InterpolationOptions
  {
  public:
    ElementLocator(const ParaFIELD& sourceField, const ProcessorGroup& distant_group);

    virtual ~ElementLocator();
    void exchangeMesh(int idistantrank,
                      MEDCouplingPointSet*& target_mesh,
                      int*& distant_ids);
    void exchangeMethod(const std::string& sourceMeth, int idistantrank, std::string& targetMeth);
    const std::vector<int>& getDistantProcIds() const { return _distant_proc_ids; }
    const MPI_Comm *getCommunicator() const;
  private:
    void _computeBoundingBoxes();
    bool _intersectsBoundingBox(int irank);
    void _exchangeMesh(MEDCouplingPointSet* local_mesh, MEDCouplingPointSet*& distant_mesh,
                       int iproc_distant, const DataArrayInt* distant_ids_send,
                       int*& distant_ids_recv);
  private:
    const ParaFIELD&  _local_para_field ;
    MEDCouplingPointSet* _local_cell_mesh;
    MEDCouplingPointSet* _local_face_mesh;
    std::vector<MEDCouplingPointSet*> _distant_cell_meshes;
    std::vector<MEDCouplingPointSet*> _distant_face_meshes;
    double* _domain_bounding_boxes;
    const ProcessorGroup& _distant_group;
    const ProcessorGroup& _local_group;
    ProcessorGroup* _union_group;
    std::vector<int> _distant_proc_ids;
  };

}

#endif
