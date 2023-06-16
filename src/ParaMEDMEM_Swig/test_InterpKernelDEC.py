#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2023  CEA/DEN, EDF R&D
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
import sys, os
import unittest
import math
from mpi4py import MPI


class ParaMEDMEM_IK_DEC_Tests(unittest.TestCase):
    """ See test_StructuredCoincidentDEC_py_1() for a quick start.
    """
    def generateFullSource(self):
        """ The complete source mesh: 4 squares each divided in 2 diagonaly (so 8 cells in total) """
        msh  = self.generateFullTarget()
        msh.simplexize(0)
        msh.setName("src_mesh")
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        fld.setMesh(msh); fld.setName("source_F");
        da = DataArrayDouble(msh.getNumberOfCells())
        da.iota()
        da *= 2
        fld.setArray(da)
        return msh, fld

    def generateFullTarget(self):
        """ The complete target mesh: 4 squares """
        m1 = MEDCouplingCMesh("tgt_msh")
        da = DataArrayDouble([0,1,2])
        m1.setCoords(da, da)
        msh = m1.buildUnstructured()
        return msh

    #
    # Below, the two functions emulating the set up of a piece of the source and target mesh
    # on each proc. Obviously in real world problems, this comes from your code and is certainly
    # not computed by cuting again from scratch the full-size mesh!!
    #
    def getPartialSource(self, rank):
        """ Will return an empty mesh piece for rank=2 and 3 """
        msh, f = self.generateFullSource()
        if rank == 0:
            sub_m, sub_f = msh[0:4], f[0:4]
        elif rank == 1:
            sub_m, sub_f = msh[4:8], f[4:8]
        sub_m.zipCoords()
        return sub_m, sub_f

    def getPartialTarget(self, rank):
        """ One square for each rank """
        msh = self.generateFullTarget()
        if rank == 2:
            sub_m = msh[[0,2]]
        elif rank == 3:
            sub_m = msh[[1,3]]
        sub_m.zipCoords()
        # Receiving side must prepare an empty field that will be filled by DEC:
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        da = DataArrayDouble(sub_m.getNumberOfCells())
        fld.setArray(da)
        fld.setName("tgt_F")
        fld.setMesh(sub_m)
        return sub_m, fld

    def testInterpKernelDEC_ctor(self):
        """ Test the various Python ctors """
        size = MPI.COMM_WORLD.size
        if size != 4:
            print("Should be run on 4 procs!")
            return
        # Define two processor groups
        nproc_source = 2
        l1, l2 = range(nproc_source), range(size - nproc_source, size)
        # With 2 iterables:
        i1 = InterpKernelDEC.New(l1, l2)
        # Should also work directly:
        i2 = InterpKernelDEC(l1, l2)
        # With 2 proc groups:
        interface = CommInterface()
        source_group = MPIProcessorGroup(interface, list(l1))
        target_group = MPIProcessorGroup(interface, list(l2))
        i3 = InterpKernelDEC.New(source_group, target_group)
        # Should also work directly:
        i4 = InterpKernelDEC(source_group, target_group)
        # With 2 iterables and a custom comm:
        i5 = InterpKernelDEC.New(l1, l2, MPI.COMM_WORLD)
        # Work directly with the **hack**
        i6 = InterpKernelDEC(l1, l2, MPI._addressof(MPI.COMM_WORLD))
        # Should fail with 2 proc groups **and** a communicator
        self.assertRaises(InterpKernelException, InterpKernelDEC.New, source_group, target_group, MPI.COMM_WORLD)
        self.assertRaises(Exception, InterpKernelDEC, source_group, target_group, MPI.COMM_WORLD)
        i6.release(); i5.release(); i4.release(); i3.release(); i2.release(); i1.release()
        source_group.release()
        target_group.release()

    @WriteInTmpDir
    def testInterpKernelDEC_2D_py_1(self):
        """ This test illustrates a basic use of the InterpKernelDEC.
        Look at the C++ documentation of the class for more informations.
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
        idec = InterpKernelDEC(source_group, target_group)

        # Write out full size meshes/fields for inspection
        if rank == 0:
            _, fld = self.generateFullSource()
            mshT = self.generateFullTarget()
            WriteField("./source_field_FULL.med", fld, True)
            WriteUMesh("./target_mesh_FULL.med", mshT, True)

        MPI.COMM_WORLD.Barrier()  # really necessary??

        #
        # OK, let's go DEC !!
        #
        if source_group.containsMyRank():
            _, fieldS = self.getPartialSource(rank)
            fieldS.setNature(IntensiveMaximum)   # The only policy supported for now ...
            WriteField("./source_field_part_%d.med" % rank, fieldS, True)
            idec.attachLocalField(fieldS)
            idec.synchronize()
            idec.sendData()

        if target_group.containsMyRank():
            mshT, fieldT = self.getPartialTarget(rank)
            fieldT.setNature(IntensiveMaximum)
            WriteUMesh("./target_mesh_part_%d.med" % rank, mshT, True)
            idec.attachLocalField(fieldT)
            idec.synchronize()
            idec.recvData()
            # Now the actual checks:
            if rank == 2:
                self.assertEqual(fieldT.getArray().getValues(), [1.0, 9.0])
            elif rank == 3:
                self.assertEqual(fieldT.getArray().getValues(), [5.0, 13.0])

        # Release DEC (this involves MPI exchanges -- notably the release of the communicator -- so better be done before MPI.Finalize()
        idec.release()
        source_group.release()
        target_group.release()
        MPI.COMM_WORLD.Barrier()

    @WriteInTmpDir
    def test_InterpKernelDEC_2D_py_2(self):
        """ More involved test using Para* objects.
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 5:
            print("Should be run on 5 procs!")
            return

        print(rank)
        nproc_source = 3
        procs_source = list(range(nproc_source))
        procs_target = list(range(size - nproc_source + 1, size))

        interface = CommInterface()
        target_group = MPIProcessorGroup(interface, procs_target)
        source_group = MPIProcessorGroup(interface, procs_source)
        dec = InterpKernelDEC(source_group, target_group)

        data_dir = os.path.join(os.environ['MEDCOUPLING_ROOT_DIR'], "share", "resources", "med")
        if not os.path.isdir(data_dir):
            data_dir = os.environ.get('MED_RESOURCES_DIR',"::").split(":")[1]

        filename_xml1 = os.path.join(data_dir, "square1_split")
        filename_xml2 = os.path.join(data_dir, "square2_split")

        MPI.COMM_WORLD.Barrier()
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
            icocofield = ICoCoMEDDoubleField(parafield.getField())
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
            icocofield = ICoCoMEDDoubleField(parafield.getField())
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

        # Some clean up that still needs MPI communication, so to be done before MPI_Finalize()
        parafield.release()
        paramesh.release()
        dec.release()
        target_group.release()
        source_group.release()
        MPI.COMM_WORLD.Barrier()

if __name__ == "__main__":
    unittest.main()
    MPI.Finalize()

