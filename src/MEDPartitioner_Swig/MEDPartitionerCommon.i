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

%module MEDPartitioner

%include std_string.i

%{
#include "MEDFileData.hxx"

#include "MEDPARTITIONER_MEDPartitioner.hxx"
#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_Graph.hxx"

using namespace MEDCoupling;
using namespace INTERP_KERNEL;
using namespace MEDPARTITIONER;
%}

%feature("autodoc", "1");
%feature("docstring");

%newobject MEDPARTITIONER::MEDPartitioner::New;
%newobject MEDPARTITIONER::MEDPartitioner::Graph;
%newobject MEDPARTITIONER::MEDPartitioner::Graph::getGraph;
%newobject MEDPARTITIONER::MEDPartitioner::Graph::getPartition;
%newobject MEDPARTITIONER::MEDPartitioner::getMEDFileData;
%feature("unref") MEDCoupling::MEDFileData "$this->decrRef();"

%nodefaultctor;

%rename (InterpKernelException) INTERP_KERNEL::Exception;

namespace MEDPARTITIONER
{
  class Graph
  {
  public:
    typedef enum {METIS,SCOTCH} splitter_type;
  public:
    virtual void partGraph(int ndomain, const std::string& options_string="", ParaDomainSelector *sel=0) throw(INTERP_KERNEL::Exception);
    MEDCoupling::MEDCouplingSkyLineArray *getGraph() const
    {
      const MEDCoupling::MEDCouplingSkyLineArray *ret(self->getGraph());
      if(ret)
        ret->incrRef();
      return const_cast<MEDCoupling::MEDCouplingSkyLineArray *>(ret);
    }
    const MEDCoupling::MEDCouplingSkyLineArray *getPartition() const
    {
      const MEDCoupling::MEDCouplingSkyLineArray *ret(self->getPartition());
      if(ret)
        ret->incrRef();
      return const_cast<MEDCoupling::MEDCouplingSkyLineArray *>(ret);
    }
    int nbVertices() const;
  };

  class MEDPartitioner
  {
  public:
    MEDPartitioner(const std::string& filename, int ndomains=1, const std::string& library="metis",bool create_boundary_faces=false, bool create_joints=false, bool mesure_memory=false) throw(INTERP_KERNEL::Exception);
    MEDPartitioner(const MEDCoupling::MEDFileData* fileData, int ndomains=1, const std::string& library="metis",bool create_boundary_faces=false, bool create_joints=false, bool mesure_memory=false) throw(INTERP_KERNEL::Exception);
    MEDPartitioner(const MEDCoupling::MEDFileData* fileData, Graph* graph, bool create_boundary_faces=false, bool create_joints=false, bool mesure_memory=false) throw(INTERP_KERNEL::Exception);
    static MEDPARTITIONER::Graph* Graph(MEDCoupling::MEDCouplingSkyLineArray* graph, Graph::splitter_type split=Graph::METIS, int* edgeweight=0) throw(INTERP_KERNEL::Exception);
    MEDCoupling::MEDFileData* getMEDFileData() throw(INTERP_KERNEL::Exception);
    void write(const std::string& filename) throw(INTERP_KERNEL::Exception);
  };
}

