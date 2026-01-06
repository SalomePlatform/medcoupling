#  -*- coding: utf-8 -*-
# Copyright (C) 2025-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

# And here we use mpi4py ability to provide its internal (C++) pointer to the communicator:
# NB: doing a proper typemap from MPI_Comm from Python to C++ requires the inclusion of mpi4py headers and .i file ... an extra dependency ...
def _IKDEC_WithComm_internal(src_procs, tgt_procs, mpicomm=None):
    from mpi4py import MPI
    import medcoupling as pmm

    # Check iterable:
    try:
        s, t = [el for el in src_procs], [el for el in tgt_procs]
    except:
        s, t = None, None
    msg = "InterpKernelDEC: invalid type in ctor arguments! Possible signatures are:\n"
    msg += "   - InterpKernelDEC(ProcessorGroup, ProcessorGroup)\n"
    msg += "   - InterpKernelDEC(<iterable>, <iterable>)\n"
    msg += "   - InterpKernelDEC(<iterable>, <iterable>, MPI_Comm*) : WARNING here the address of the communicator should be passed with MPI._addressof(the_com)\n"
    msg += "   - InterpKernelDEC.New(ProcessorGroup, ProcessorGroup)\n"
    msg += "   - InterpKernelDEC.New(<iterable>, <iterable>)\n"
    msg += "   - InterpKernelDEC.New(<iterable>, <iterable>, MPI_Comm)\n"
    if mpicomm is None:
        if isinstance(src_procs, pmm.ProcessorGroup) and isinstance(
            tgt_procs, pmm.ProcessorGroup
        ):
            return pmm.InterpKernelDEC._NewWithPG_internal(src_procs, tgt_procs)
        elif not s is None:  # iterable
            return pmm.InterpKernelDEC._NewWithComm_internal(
                s, t, MPI._addressof(MPI.COMM_WORLD)
            )
        else:
            raise pmm.InterpKernelException(msg)
    else:
        if s is None:
            raise pmm.InterpKernelException(msg)  # must be iterable
        return pmm.InterpKernelDEC._NewWithComm_internal(s, t, MPI._addressof(mpicomm))


def _ODEC_WithComm_internal(procs, mpicomm=None):
    from mpi4py import MPI
    import medcoupling as pmm

    # Check iterable:
    try:
        g = [el for el in procs]
    except:
        msg = "OverlapDEC: invalid type in ctor arguments! Possible signatures are:\n"
        msg += "   - OverlapDEC.New(<iterable>)\n"
        msg += "   - OverlapDEC.New(<iterable>, MPI_Comm)\n"
        msg += "   - OverlapDEC(<iterable>)\n"
        msg += "   - OverlapDEC(<iterable>, MPI_Comm*) : WARNING here the address of the communicator should be passed with MPI._addressof(the_com)\n"
        raise pmm.InterpKernelException(msg)
    if mpicomm is None:
        return pmm.OverlapDEC(g)
    else:
        return pmm.OverlapDEC._NewWithComm_internal(g, MPI._addressof(mpicomm))
