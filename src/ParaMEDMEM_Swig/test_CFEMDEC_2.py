#!/usr/bin/env python3
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

"""
See EDF35188 : Test P1P1 on 3 procs. 1 source procs and1 target procs.

This test check the correction of the remapper FE
"""

import medcoupling as mc
from mpi4py import MPI

globalComm = MPI.COMM_WORLD

size = globalComm.size
rank = globalComm.rank

if size != 3:
    raise RuntimeError("Expected to be lanched with 3 procs !")

procs_source = [
    0,
    1,
]
procs_target = [
    2,
]

idec = mc.CFEMDEC(procs_source, procs_target)

mesh_name = "00000003"
field_name = "00000008"

if rank in procs_source:
    file_name = f"test_CFEMDEC_2_pres_{rank}.med"

    src_mesh = mc.MEDFileMesh.New(file_name, mesh_name)
    src_mesh_interf = src_mesh[0]
    src_globalNodeIds = src_mesh.getGlobalNumFieldAtLevel(1)
    src_field = mc.ReadFieldNode(file_name, mesh_name, 0, field_name, -1, -1)
    # maybe better way to do it
    src_field.setMesh(src_mesh_interf)

    assert src_globalNodeIds.getNbOfElems() == [166, 158][rank]

    idec.attachLocalMesh(src_mesh_interf, src_globalNodeIds)
    idec.sendToTarget(src_field)

if rank in procs_target:

    file_name = "test_CFEMDEC_2_pres_SEQ.med"
    trg_mesh = mc.MEDFileMesh.New(file_name, mesh_name)
    trg_mesh_interf = trg_mesh[0]
    trg_globalNodeIds = trg_mesh.getGlobalNumFieldAtLevel(1)
    trg_field = mc.ReadFieldNode(file_name, mesh_name, 0, field_name, -1, -1)

    assert trg_globalNodeIds.getNbOfElems() == 280

    idec.attachLocalMesh(trg_mesh_interf, trg_globalNodeIds)
    trg_recv = idec.receiveFromSource()

    # compare values with reference
    assert trg_recv.getArray().isEqualWithoutConsideringStr(trg_field.getArray(), 1e-11)

# mpirun -np 3 xterm -e "gdb -x cmd.gdb --args python3 test_CFEMDEC_2.py"
