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

#ifndef __MEDPARTITIONER_MEDPARTITIONER_HXX__
#define __MEDPARTITIONER_MEDPARTITIONER_HXX__

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_Graph.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class MEDFileData;
}

namespace MEDPARTITIONER
{
  class Topology;
  class MeshCollection;

  class MEDPARTITIONER_EXPORT MEDPartitioner
  {
  public:
    MEDPartitioner(const std::string& filename, int ndomains=1, const std::string& library="metis",bool create_boundary_faces=false, bool create_joints=false, bool mesure_memory=false);
    MEDPartitioner(const MEDCoupling::MEDFileData* fileData, int ndomains=1, const std::string& library="metis",bool create_boundary_faces=false, bool create_joints=false, bool mesure_memory=false);
    MEDPartitioner(const MEDCoupling::MEDFileData* fileData, Graph* graph, bool create_boundary_faces=false, bool create_joints=false, bool mesure_memory=false);
    static MEDPARTITIONER::Graph* Graph(MEDCoupling::MEDCouplingSkyLineArray* graph, Graph::splitter_type split=Graph::METIS, int* edgeweight=0);
    void write(const std::string& filename);
    MEDCoupling::MEDFileData* getMEDFileData();
    ~MEDPartitioner();

    MEDCoupling::MEDFileData *convertToMEDFileData(MeshCollection* meshcollection);
    void createPartitionCollection(int ndomains, const std::string& library,bool create_boundary_faces, bool create_joints, bool mesure_memory);

  private:
    MeshCollection* _input_collection;
    MeshCollection* _output_collection;
    Topology*       _new_topology;
  };
}
#endif
