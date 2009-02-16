#!/usr/bin/env python

# Copyright (C) 2005  OPEN CASCADE, CEA, EDF R&D, LEG
#           PRINCIPIA R&D, EADS CCR, Lip6, BV, CEDRAT
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either 
# version 2.1 of the License.
# 
# This library is distributed in the hope that it will be useful 
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

# testIntersectionDEC_2D

MPI_Init(sys.argv)

size = MPI_Comm_size(MPI_COMM_WORLD)
rank = MPI_Comm_rank(MPI_COMM_WORLD)
if size != 5:
    raise RuntimeError, "Expect MPI_COMM_WORLD size == 5"

nproc_source = 3
procs_source = range( nproc_source )
procs_target = range( size - nproc_source + 1, size)

interface = CommInterface()

target_group = MPIProcessorGroup(interface, procs_target)
source_group = MPIProcessorGroup(interface, procs_source)

dec = IntersectionDEC(source_group, target_group)

mesh       =0
support    =0
paramesh   =0
parafield  =0
parasupport=0
icocofield =0

data_dir = os.environ['MED_ROOT_DIR']
tmp_dir  = os.environ['TMP']

if not tmp_dir or len(tmp_dir)==0:
    tmp_dir = "/tmp"
    pass

filename_xml1 = os.path.join(data_dir, "share/salome/resources/med/square1_split")
filename_xml2 = os.path.join(data_dir, "share/salome/resources/med/square2_split")

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

if source_group.containsMyRank():
    field_before_int = parafield.getVolumeIntegral(1)
    dec.synchronize()
    print "DEC usage"
    dec.setForcedRenormalization(False)

    dec.sendData()
    print "writing 1"
    paramesh.write(MED_DRIVER,"./sourcesquareb");
    filename = "./sourcesquareb_" + str(source_group.myRank()+1)
    parafield.write(MED_DRIVER, "./sourcesquareb", "boundary");

    print "b dec.recvData()"
    dec.recvData()
    print "writing 2"
    paramesh.write(MED_DRIVER, "./sourcesquare")
    parafield.write(MED_DRIVER, "./sourcesquare", "boundary")

    filename = "./sourcesquare_" + str(source_group.myRank()+1)
    field_after_int = parafield.getVolumeIntegral(1)
    print "field_before_int", field_before_int,"field_after_int", field_after_int

    pass

if target_group.containsMyRank():
    dec.synchronize()
    print "TARGET: after dec.synchronize()"
    dec.setForcedRenormalization(False)

    dec.recvData()
    paramesh.write(0, "./targetsquareb")
    parafield.write(0, "./targetsquareb", "boundary")
    dec.sendData();
    paramesh.write(0, "./targetsquare")
    parafield.write(0, "./targetsquare", "boundary")
    pass


MPI_Barrier(MPI_COMM_WORLD)
MPI_Finalize()
print "### End of testIntersectionDEC_2D"
