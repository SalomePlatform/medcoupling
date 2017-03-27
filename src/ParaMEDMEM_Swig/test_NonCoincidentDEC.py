#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

from ParaMEDMEM import *
import sys, os

MPI_Init(sys.argv)

size = MPI_Comm_size(MPI_COMM_WORLD)
rank = MPI_Comm_rank(MPI_COMM_WORLD)
if size != 5:
    raise RuntimeError("Expect MPI_COMM_WORLD size == 5")

nproc_source = 3
procs_source = list(range(nproc_source))
procs_target = list(range(size - nproc_source + 1, size))

interface = CommInterface()

target_group = MPIProcessorGroup(interface, procs_target)
source_group = MPIProcessorGroup(interface, procs_source)

source_mesh= 0
target_mesh= 0
parasupport= 0
mesh       = 0
support    = 0
field      = 0
paramesh   = 0
parafield  = 0
icocofield = 0

dec = NonCoincidentDEC(source_group, target_group)

data_dir = os.environ['MEDCOUPLING_ROOT_DIR']
tmp_dir  = os.environ['TMP']
if tmp_dir == '':
    tmp_dir = "/tmp"
    pass

filename_xml1 = data_dir + "/share/resources/med/square1_split"
filename_xml2 = data_dir + "/share/resources/med/square2_split"

MPI_Barrier(MPI_COMM_WORLD)

if source_group.containsMyRank():

    filename = filename_xml1 + str(rank+1) + ".med"
    meshname = "Mesh_2_" + str(rank+1)

    mesh = MESH(MED_DRIVER, filename, meshname)
    support = SUPPORT(mesh, "all elements", MED_CELL)
    paramesh = ParaMESH(mesh, source_group, "source mesh")

    parasupport = UnstructuredParaSUPPORT( support, source_group)
    comptopo = ComponentTopology()

    parafield = ParaFIELD(parasupport, comptopo)

    nb_local = support.getNumberOfElements(MED_ALL_ELEMENTS);

    value = [1.0]*nb_local

    parafield.getField().setValue(value)
    icocofield = ICoCo_MEDField(paramesh,parafield)
    dec.attachLocalField(icocofield,'P0')
    pass

if target_group.containsMyRank():

    filename = filename_xml2 + str(rank - nproc_source + 1) + ".med"
    meshname = "Mesh_3_" + str(rank - nproc_source + 1)

    mesh = MESH(MED_DRIVER, filename, meshname)
    support = SUPPORT(mesh, "all elements", MED_CELL)
    paramesh = ParaMESH(mesh, target_group, "target mesh")

    parasupport = UnstructuredParaSUPPORT( support, target_group)
    comptopo = ComponentTopology()
    parafield = ParaFIELD(parasupport, comptopo)

    nb_local = support.getNumberOfElements(MED_ALL_ELEMENTS)
    value = [0.0]*nb_local

    parafield.getField().setValue(value)
    icocofield = ICoCo_MEDField(paramesh,parafield)

    dec.attachLocalField(icocofield, 'P0')
    pass

field_before_int = [0.0]
field_after_int = [0.0]

if source_group.containsMyRank():

    field_before_int = [parafield.getVolumeIntegral(1)]
    MPI_Bcast(field_before_int, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    dec.synchronize()
    print("DEC usage")
    dec.setForcedRenormalization(False)

    dec.sendData()
    pass

if target_group.containsMyRank():

    MPI_Bcast(field_before_int, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD)
    dec.synchronize()
    dec.setForcedRenormalization(False)
    dec.recvData()
    field_after_int = [parafield.getVolumeIntegral(1)]
    pass

MPI_Bcast(field_before_int, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD)
MPI_Bcast(field_after_int , 1, MPI_DOUBLE, size-1, MPI_COMM_WORLD)

epsilon = 1e-6
if abs(field_before_int[0] - field_after_int[0]) > epsilon:
    print("Field before is not equal field after: %s != %s"%\
          (field_before_int[0],field_after_int[0]))
    pass


MPI_Barrier(MPI_COMM_WORLD)
MPI_Finalize()
print("# End of testNonCoincidentDEC")
