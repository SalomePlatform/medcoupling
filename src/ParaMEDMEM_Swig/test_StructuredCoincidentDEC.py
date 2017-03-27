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

class ParaMEDMEMBasicsTest2(unittest.TestCase):
    def testStructuredCoincidentDEC(self):
        MPI_Init(sys.argv)
        #
        size = MPI_Comm_size(MPI_COMM_WORLD)
        rank = MPI_Comm_rank(MPI_COMM_WORLD)
        #
        if size < 4:
            raise RuntimeError("Expect MPI_COMM_WORLD size >= 4")
        #
        interface = CommInterface()
        #
        self_group   = MPIProcessorGroup(interface, rank, rank)
        target_group = MPIProcessorGroup(interface, 3, size-1)
        source_group = MPIProcessorGroup(interface, 0, 2)
        #
        mesh      = 0
        support   = 0
        paramesh  = 0
        parafield = 0
        comptopo  = 0
        icocofield= 0
        #
        data_dir = os.environ['MEDCOUPLING_ROOT_DIR']
        tmp_dir  = os.environ['TMP']
        if tmp_dir == '':
            tmp_dir = "/tmp"
            pass

        filename_xml1    = data_dir + "/share/resources/med/square1_split"
        filename_2       = data_dir + "/share/resources/med/square1.med"
        filename_seq_wr  = tmp_dir + "/"
        filename_seq_med = tmp_dir + "/myWrField_seq_pointe221.med"

        dec = StructuredCoincidentDEC(source_group, target_group)
        MPI_Barrier(MPI_COMM_WORLD)
        if source_group.containsMyRank():
            filename = filename_xml1 + str(rank+1) + ".med"
            meshname = "Mesh_2_" + str(rank+1)
            mesh=ReadUMeshFromFile(filename,meshname,0)
            paramesh=ParaMESH(mesh,source_group,"source mesh")
            comptopo=ComponentTopology(6)
            parafield=ParaFIELD(ON_CELLS,NO_TIME,paramesh,comptopo)
            parafield.getField().setNature(IntensiveMaximum)
            nb_local=mesh.getNumberOfCells()
            global_numbering=paramesh.getGlobalNumberingCell2()
            value = []
            for ielem in range(nb_local):
                for icomp in range(6):
                    value.append(global_numbering[ielem]*6.0+icomp);
                    pass
                pass
            parafield.getField().setValues(value)
            icocofield = ICoCoMEDField(mesh,parafield.getField())
            dec.setMethod("P0")
            dec.attachLocalField(parafield)
            dec.synchronize()
            dec.sendData()
            pass

        if target_group.containsMyRank():
            meshname2 = "Mesh_2"
            mesh=ReadUMeshFromFile(filename_2, meshname2,0)
            paramesh=ParaMESH(mesh, self_group, "target mesh")
            comptopo=ComponentTopology(6,target_group)
            parafield=ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo)
            parafield.getField().setNature(IntensiveMaximum)
            nb_local=mesh.getNumberOfCells()
            value = [0.0]*(nb_local*comptopo.nbLocalComponents())
            parafield.getField().setValues(value)
            icocofield = ICoCoMEDField(mesh,parafield.getField())
            dec.setMethod("P0")
            dec.attachLocalField(parafield)
            dec.synchronize()
            dec.recvData()
            recv_value = parafield.getField().getArray().getValues()
            for i in range(nb_local):
                first=comptopo.firstLocalComponent()
                for icomp in range(comptopo.nbLocalComponents()):
                    self.assertTrue(math.fabs(recv_value[i*comptopo.nbLocalComponents()+icomp]-
                                              (float)(i*6+icomp+first))<1e-12)
                    pass
                pass
            pass
        comptopo=0
        interface = 0
        mesh       =0
        support    =0
        paramesh   =0
        parafield  =0
        icocofield =0
        dec=0
        self_group =0
        target_group = 0
        source_group = 0
        MPI_Barrier(MPI_COMM_WORLD)
        MPI_Finalize()
        print("End of test StructuredCoincidentDEC")
        pass


unittest.main()
