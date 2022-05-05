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

#ifndef __INTERPKERNELDEC_HXX__
#define __INTERPKERNELDEC_HXX__

#include "DisjointDEC.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationOptions.hxx"

namespace MEDCoupling
{
  class InterpolationMatrix;

  /*!
    \anchor InterpKernelDEC-det
    \class InterpKernelDEC

    \section InterpKernelDEC-over Overview

    The InterpKernelDEC enables the \ref InterpKerRemapGlobal "remapping" (or interpolation) of fields between
    two parallel codes.

    The projection
    methodology is based on the algorithms of %INTERP_KERNEL, that is to say, they work in a similar fashion than
    what the \ref remapper "sequential remapper" does. The following \ref discretization "projection methods"
    are supported: P0->P0 (the most common case), P1->P0, P0->P1.

    The computation is possible for 3D meshes, 2D meshes, and 3D-surface
    meshes. Dimensions must be identical for code A and code B (for instance, though it could be
    desirable, it is not yet possible to couple 3D surfaces with 2D surfaces).

    The name "InterpKernelDEC" comes from the fact that this class uses exactly the same algorithms
    as the sequential remapper. Both this class and the sequential
    \ref MEDCoupling::MEDCouplingRemapper "MEDCouplingRemapper" are built on top of the %INTERP_KERNEL
    algorithms (notably the computation of the intersection volumes).

    Among the important properties inherited from the parent abstract class \ref DisjointDEC-det "DisjointDEC",
    the two \ref MPIProcessorGroup-det "processor groups" (source and target) must have a void intersection.

    \image html NonCoincident_small.png "Transfer of a field supported by a quadrangular mesh to a triangular mesh".

    \image latex NonCoincident_small.eps "Transfer of a field supported by a quadrangular mesh to a triangular mesh"

    In the figure above we see the transfer of a field based on a quadrangular mesh to a new field supported by
    a triangular mesh. In a P0-P0 interpolation, to obtain the value on a triangle, the values on the
    quadrangles are weighted by their intersection area and summed.

    A typical use of InterpKernelDEC encompasses two distinct phases :
    - A setup phase during which the intersection volumes are computed and the communication structures are
    setup. This corresponds to calling the InterpKernelDEC::synchronize() method.
    - A running phase during which the projections are actually performed. This corresponds to the calls to
    sendData() and recvData() which actually trigger the data exchange. The data exchange are synchronous
    in the current version of the library so that recvData() and sendData() calls must be synchronized
    on code A and code B processor groups.

    The following code excerpt illustrates a typical use of the InterpKernelDEC class.

    \code
    ...
    InterpKernelDEC dec(groupA, groupB);
    dec.attachLocalField(field);
    dec.synchronize();
    if (groupA.containsMyRank())
    dec.recvData();
    else if (groupB.containsMyRank())
    dec.sendData();
    ...
    \endcode
    A \ref InterpKerRemapGlobal "remapping" of the field from the source mesh to the target mesh is performed by
    the function synchronise(), which computes the interpolation matrix.

    Computing the field on the receiving side can be expressed in terms of a matrix-vector product :
    \f$ \phi_t=W.\phi_s\f$, with \f$ \phi_t \f$ the field on the target side and \f$ \phi_s \f$ the field
    on the source side.
    When remapping a 3D surface to another 3D surface, a projection phase is necessary to match elements
    from both sides. Care must be taken when defining this projection to obtain a
    \ref InterpKerRemapGlobal "conservative remapping".

    In the P0-P0 case, this matrix is a plain rectangular matrix with coefficients equal to the
    intersection areas between triangle and quadrangles. For instance, in the above figure, the matrix
    is :

    \f[
    \begin{tabular}{|cccc|}
    0.72 & 0 & 0.2 & 0 \\
    0.46 & 0 & 0.51 & 0.03\\
    0.42 & 0.53 & 0 & 0.05\\
    0 & 0 & 0.92 & 0.05 \\
    \end{tabular}
    \f]

    \section InterpKernelDEC-options Options
    On top of the usual \ref MEDCoupling::DECOptions "DEC options", the options supported by %InterpKernelDEC objects are
    related to the underlying \ref InterpKerIntersectors "intersector class".
    All the options available in the intersector objects are
    available for the %InterpKernelDEC object. The various options available for  intersectors can
    be reviewed in \ref InterpKerIntersectors.

    For instance :
    \verbatim
    InterpKernelDEC dec(source_group, target_group);
    dec.attachLocalField(field);
    dec.setDoRotate(false);
    dec.setPrecision(1e-12);
    dec.synchronize();
    \endverbatim

    \warning{  Options must be set before calling the synchronize method. }
  */

  class InterpKernelDEC : public DisjointDEC, public INTERP_KERNEL::InterpolationOptions
  {
  public:  
    InterpKernelDEC();
    InterpKernelDEC(ProcessorGroup& source_group, ProcessorGroup& target_group);
    InterpKernelDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids,
                    const MPI_Comm& world_comm=MPI_COMM_WORLD);
    virtual ~InterpKernelDEC();
    void release();

    void synchronize();
    void recvData();
    void recvData(double time);
    void sendData();
    void sendData(double time , double deltatime);
    void prepareSourceDE() { }
    void prepareTargetDE() { }
  private :
    InterpolationMatrix* _interpolation_matrix;
  };
}

#endif
