// Copyright (C) 2007-2023  CEA/DEN, EDF R&D
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

#ifndef __STRUCTUREDCOINCIDENTDEC_HXX__
#define __STRUCTUREDCOINCIDENTDEC_HXX__

#include "DisjointDEC.hxx"
#include "BlockTopology.hxx"


namespace MEDCoupling
{
  class DEC;
  class BlockTopology;

  /*!
    \anchor StructuredCoincidentDEC-det
    \class StructuredCoincidentDEC

    This class aims at \ref interpolation "remapping fields" that have identical
    structured supports (=the same underlying mesh) but different parallel topologies
    (=different sub-domains in the structured mesh). It can be used to couple
    together multi-physics codes that operate on the same domain
    with different partitioning. This can be useful for example if one of
    the computation is much faster than the other. It can also be used
    to couple together codes that share an interface that was generated
    in the same manner (with identical global ids).
    Also, this \ref para-dec "DEC" can be used for fields that have component topologies,
    i.e., components that are scattered over several processors.

    The remapping between the two supports is based on identity of global
    ids, instead of geometrical considerations (as it is the case for
    \ref InterpKernelDEC-det "InterpKernelDEC").
    Therefore, beware that this \ref para-dec "DEC" can not be used
    for coincident meshes if they do *not* have the exact same numbering.

    With this %DEC no projection, and no interpolation of the field data is done, contrary
    to what happens in \ref InterpKernelDEC-det "InterpKernelDEC". It is just
    a matter of allocating the values from one side to the other, using directly the cell
    identifiers.

    As all the other DECs, its usage requires two phases :
    - a setup phase during which the topologies are exchanged so that
    the target side knows from which processors it should expect
    the data.
    - a send/recv phase during which the field data is actually transferred.

    This example illustrates the sending of a field with
    the \c StructuredCoincidentDEC :
    \code
    ...
    StructuredCoincidentDEC dec(groupA, groupB);
    dec.attachLocalField(field);
    dec.synchronize();
    if (groupA.containsMyRank())
    dec.recvData();
    else if (groupB.containsMyRank())
    dec.sendData();
    ...
    \endcode

    Creating a ParaFIELD to be attached to the %DEC is done in exactly the same way as for
    the other DECs, if only the partitioning of the support mesh differs.
    In the case where the
    fields have also different *component* topologies, creating the ParaFIELD
    requires some more effort. See the \ref para-over "parallelism" section for more details.
  */
  class StructuredCoincidentDEC : public DisjointDEC
  {
  public:
    StructuredCoincidentDEC();
    StructuredCoincidentDEC( ProcessorGroup& source, ProcessorGroup& target);
    virtual ~StructuredCoincidentDEC();
    void release();

    void synchronize();
    void recvData();
    void sendData();
    void prepareSourceDE();
    void prepareTargetDE();

  private :
    void synchronizeTopology();
    void broadcastTopology(BlockTopology*&, int tag);

    BlockTopology* _topo_source;
    BlockTopology* _topo_target;

    bool _owns_topo_source;
    bool _owns_topo_target;

    int* _send_counts;
    int* _recv_counts;
    int* _send_displs;
    int* _recv_displs;
    double* _recv_buffer;
    double* _send_buffer;
  };
}

#endif
