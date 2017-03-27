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
import unittest
import math

class ParaMEDMEMBasicsTest(unittest.TestCase):
    def testInterpKernelDEC_2D(self):
        MPI_Init(sys.argv)
        size = MPI_Comm_size(MPI_COMM_WORLD)
        rank = MPI_Comm_rank(MPI_COMM_WORLD)
        if size != 5:
            raise RuntimeError("Expect MPI_COMM_WORLD size == 5")
        print(rank)
        nproc_source = 3
        procs_source = list(range(nproc_source))
        procs_target = list(range(size - nproc_source + 1, size))

        interface = CommInterface()
        target_group = MPIProcessorGroup(interface, procs_target)
        source_group = MPIProcessorGroup(interface, procs_source)
        dec = InterpKernelDEC(source_group, target_group)

        mesh       =0
        support    =0
        paramesh   =0
        parafield  =0
        icocofield =0
        data_dir = os.environ['MEDCOUPLING_ROOT_DIR']
        tmp_dir  = os.environ['TMP']

        if not tmp_dir or len(tmp_dir)==0:
            tmp_dir = "/tmp"
            pass

        filename_xml1 = os.path.join(data_dir, "share/resources/med/square1_split")
        filename_xml2 = os.path.join(data_dir, "share/resources/med/square2_split")

        MPI_Barrier(MPI_COMM_WORLD)
        if source_group.containsMyRank():
            filename = filename_xml1 + str(rank+1) + ".med"
            meshname = "Mesh_2_" + str(rank+1)
            mesh=ReadUMeshFromFile(filename,meshname,0)
            paramesh=ParaMESH(mesh,source_group,"source mesh")
            comptopo = ComponentTopology()
            parafield = ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo)
            parafield.getField().setNature(IntensiveMaximum)
            nb_local=mesh.getNumberOfCells()
            value = [1.0]*nb_local
            parafield.getField().setValues(value)
            icocofield = ICoCoMEDField(mesh,parafield.getField())
            dec.setMethod("P0")
            dec.attachLocalField(icocofield)
            pass
        else:
            filename = filename_xml2 + str(rank - nproc_source + 1) + ".med"
            meshname = "Mesh_3_" + str(rank - nproc_source + 1)
            mesh=ReadUMeshFromFile(filename,meshname,0)
            paramesh=ParaMESH(mesh,target_group,"target mesh")
            comptopo = ComponentTopology()
            parafield = ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo)
            parafield.getField().setNature(IntensiveMaximum)
            nb_local=mesh.getNumberOfCells()
            value = [0.0]*nb_local
            parafield.getField().setValues(value)
            icocofield = ICoCoMEDField(mesh,parafield.getField())
            dec.setMethod("P0")
            dec.attachLocalField(icocofield)
            pass

        if source_group.containsMyRank():
            field_before_int = parafield.getVolumeIntegral(0,True)
            dec.synchronize()
            dec.setForcedRenormalization(False)
            dec.sendData()
            dec.recvData()
            field_after_int=parafield.getVolumeIntegral(0,True);
            self.assertTrue(math.fabs(field_after_int-field_before_int)<1e-8)
            pass
        else:
            dec.synchronize()
            dec.setForcedRenormalization(False)
            dec.recvData()
            dec.sendData()
            pass
        ## end
        interface = 0
        target_group = 0
        source_group = 0
        dec = 0
        mesh       =0
        support    =0
        paramesh   =0
        parafield  =0
        icocofield =0
        MPI_Barrier(MPI_COMM_WORLD)
        MPI_Finalize()
        pass
    pass

unittest.main()
