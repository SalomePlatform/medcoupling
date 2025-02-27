#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2025  CEA, EDF
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

from medcoupling import *
from ParaMEDMEMTestTools import WriteInTmpDir
import os
import unittest
import math
from mpi4py import MPI


class ParaMEDMEM_SC_DEC_Tests(unittest.TestCase):
    """ See test_StructuredCoincidentDEC_py_1() for a quick start.
    """

    def generateFullMeshField(self):
        """ The complete mesh: 4 squares each divided in 2 diagonaly (so 8 cells in total)
        Note that in this case, this is the **only** mesh for the whole problem.
        """
        m1 = MEDCouplingCMesh("tgt_msh")
        da = DataArrayDouble([0,1,2])
        m1.setCoords(da, da)
        msh = m1.buildUnstructured()

        msh.simplexize(0)
        msh.setName("src_mesh")
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        fld.setMesh(msh); fld.setName("source_F")
        da = DataArrayDouble(msh.getNumberOfCells())
        da.iota()
        da *= 2
        fld.setArray(da)
        return msh, fld

    #
    # Below, the function emulating the set up of a piece of the mesh being owned by
    # a given processor.
    #
    def getPartialSource(self, rank):
        msh, f = self.generateFullMeshField()
        if rank == 0:
            sub_ids = [0,1,4,5]
        elif rank == 1:
            sub_ids = [2,3,6,7]
        sub_m, sub_f = msh[sub_ids], f[sub_ids]
        sub_m.zipCoords()
        return sub_m, sub_f

    def getPartialTarget(self, rank):
        msh, f = self.generateFullMeshField()
        if rank == 2:
            sub_ids = [0,1,2,3]
        elif rank == 3:
            sub_ids = [4,5,6,7]
        sub_m, sub_f = msh[sub_ids], f[sub_ids]
        sub_m.zipCoords()
        return sub_m, sub_f

    @WriteInTmpDir
    def test_StructuredCoincidentDEC_py_1(self):
        """ This test illustrates a basic use of the StructuredCoincidentDEC which allows to
        resdistribute a field/mesh which is already scattered on several processors into a different configuration.
        Look at the C++ documentation of the class for more informations.
        Note that in the case of the StructuredCoincidentDEC no interpolation whatsoever is performed. This is only
        really a redistribution of the data among the processors.
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 4:
            print("Should be run on 4 procs!")
            return

        # Define two processor groups
        nproc_source = 2
        procs_source = list(range(nproc_source))
        procs_target = list(range(size - nproc_source, size))

        interface = CommInterface()
        source_group = MPIProcessorGroup(interface, procs_source)
        target_group = MPIProcessorGroup(interface, procs_target)

        scdec = StructuredCoincidentDEC(source_group, target_group)

        # Write out full size meshes/fields for inspection
        if rank == 0:
            _, fld = self.generateFullMeshField()
            WriteField("./source_field_FULL.med", fld, True)

        #
        # OK, let's go DEC !!
        #
        if source_group.containsMyRank():
            _, fieldS = self.getPartialSource(rank)
            fieldS.setNature(IntensiveMaximum)   # The only policy supported for now ...
            WriteField("./source_field_part_%d.med" % rank, fieldS, True)
            scdec.attachLocalField(fieldS)
            scdec.synchronize()
            scdec.sendData()

        if target_group.containsMyRank():
            mshT, fieldT = self.getPartialTarget(rank)
            fieldT.setNature(IntensiveMaximum)
            WriteUMesh("./target_mesh_part_%d.med" % rank, mshT, True)
            scdec.attachLocalField(fieldT)
            scdec.synchronize()
            scdec.recvData()
            # Now the actual checks:
            if rank == 2:
                self.assertEqual(fieldT.getArray().getValues(), [0.0, 2.0, 8.0, 10.0])
            elif rank == 3:
                self.assertEqual(fieldT.getArray().getValues(), [4.0, 6.0, 12.0, 14.0])

        # Release DEC (this involves MPI exchanges -- so better be done before MPI.Finalize()
        scdec.release()
        source_group.release()
        target_group.release()

        MPI.COMM_WORLD.Barrier()

    @WriteInTmpDir
    def test_StructuredCoincidentDEC_py_2(self):
        """ More involved tests using Para* objects ...
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank

        if size < 4:
            print("Should be run on >= 4 procs!")
            return

        interface = CommInterface()

        source_group = MPIProcessorGroup(interface, 0, 2)
        target_group = MPIProcessorGroup(interface, 3, size-1)
        self_group = MPIProcessorGroup(interface, rank, rank)

        data_dir = os.path.join(os.environ['MEDCOUPLING_ROOT_DIR'], "share", "resources", "med")
        if not os.path.isdir(data_dir):
            data_dir = os.environ.get('MED_RESOURCES_DIR',"::").split(":")[1]
        tmp_dir  = os.environ.get('TMP', "")
        if tmp_dir == '':
            tmp_dir = "/tmp"

        filename_xml1 = os.path.join(data_dir, "square1_split")
        filename_2 = os.path.join(data_dir, "square1.med")
        filename_seq_wr  = "."
        filename_seq_med = os.path.join(".", "myWrField_seq_pointe221.med")

        dec = StructuredCoincidentDEC(source_group, target_group)

        MPI.COMM_WORLD.Barrier()
        if source_group.containsMyRank():
            filename = filename_xml1 + str(rank+1) + ".med"
            meshname = "Mesh_2_" + str(rank+1)
            mesh = ReadUMeshFromFile(filename,meshname,0)
            paramesh=ParaMESH(mesh,source_group,"source mesh")
            comptopo=ComponentTopology(6)
            parafield=ParaFIELD(ON_CELLS,NO_TIME,paramesh,comptopo)
            parafield.getField().setNature(IntensiveMaximum)
            nb_local=mesh.getNumberOfCells()
            global_numbering=paramesh.getGlobalNumberingCell2()
            value = []
            for ielem in range(nb_local):
                for icomp in range(6):
                    value.append(global_numbering[ielem]*6.0+icomp)
            parafield.getField().setValues(value)
            icocofield = ICoCoMEDDoubleField(parafield.getField())
            dec.setMethod("P0")
            dec.attachLocalField(parafield)
            dec.synchronize()
            dec.sendData()

        if target_group.containsMyRank():
            meshname2 = "Mesh_2"
            mesh = ReadUMeshFromFile(filename_2, meshname2,0)
            paramesh = ParaMESH(mesh, self_group, "target mesh")
            comptopo = ComponentTopology(6,target_group)
            parafield = ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo)
            parafield.getField().setNature(IntensiveMaximum)
            nb_local=mesh.getNumberOfCells()
            value = [0.0]*(nb_local*comptopo.nbLocalComponents())
            parafield.getField().setValues(value)
            icocofield = ICoCoMEDDoubleField(parafield.getField())
            dec.setMethod("P0")
            dec.attachLocalField(parafield)
            dec.synchronize()
            dec.recvData()
            recv_value = parafield.getField().getArray().getValues()
            for i in range(nb_local):
                first=comptopo.firstLocalComponent()
                for icomp in range(comptopo.nbLocalComponents()):
                    self.assertTrue(math.fabs(recv_value[i*comptopo.nbLocalComponents()+icomp]-(float)(i*6+icomp+first))<1e-12)

        # Release DEC (this involves MPI exchanges -- so better be done before MPI.Finalize()
        parafield.release()
        paramesh.release()
        dec.release()
        target_group.release()
        source_group.release()
        self_group.release()

        MPI.COMM_WORLD.Barrier()

if __name__ == "__main__":
    unittest.main()
    MPI.Finalize()
