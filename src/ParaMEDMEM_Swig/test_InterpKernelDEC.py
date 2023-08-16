#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2023  CEA, EDF
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

    def test_InterpKernelDEC_default(self):
        """
        [EDF27375] : Put a default value when non intersecting case
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 4:
            print("Should be run on 4 procs!")
            return
        nproc_source = 2
        procs_source = list(range(nproc_source))
        procs_target = list(range(size - nproc_source, size))

        interface = CommInterface()
        target_group = MPIProcessorGroup(interface, procs_target)
        source_group = MPIProcessorGroup(interface, procs_source)
        dec = InterpKernelDEC(source_group, target_group)
        #
        MPI.COMM_WORLD.Barrier()
        if source_group.containsMyRank():
            mesh = eval("Source_Proc_{}".format(rank))()
            nb_local=mesh.getNumberOfCells()
            field = MEDCouplingFieldDouble(ON_CELLS)
            field.setNature(IntensiveMaximum)
            field.setMesh(mesh)
            arr = DataArrayDouble(nb_local) ; arr.iota() ; arr += rank
            field.setArray(arr)
            dec.attachLocalField(field)
            dec.synchronizeWithDefaultValue(-2000.0)
            dec.sendData()
            # target -> source
            dec.recvData()
            if rank == 0:
                self.assertTrue(field.getArray().isEqual(DataArrayDouble([0.6,0.6,-2000]),1e-12))
                self.assertTrue( dec.retrieveNonFetchedIds().isEqual(DataArrayInt([2])) )
            if rank == 1:
                self.assertTrue(field.getArray().isEqual(DataArrayDouble([1.0,-2000]),1e-12))
                self.assertTrue( dec.retrieveNonFetchedIds().isEqual(DataArrayInt([1])) )
        else:
            mesh = eval("Target_Proc_{}".format(rank))()
            nb_local=mesh.getNumberOfCells()
            field = MEDCouplingFieldDouble(ON_CELLS)
            field.setNature(IntensiveMaximum)
            field.setMesh(mesh)
            arr = DataArrayDouble(nb_local) ; arr[:] = -20
            field.setArray(arr)
            dec.attachLocalField(field)
            dec.synchronizeWithDefaultValue(-1000.0)
            dec.recvData()
            if rank == 2:
                # matrix S0 / T2 = [[(0,S0,1),(1,S0,1.5)]
                # IntensiveMaximum => [[(0,S0,1/2.5),(1,S0,1.5/2.5)]
                # 
                self.assertTrue(field.getArray().isEqual(DataArrayDouble([0.6]),1e-12))
                self.assertTrue( dec.retrieveNonFetchedIds().isEqual(DataArrayInt([])) )
            if rank == 3:
                # matrix S1 / T3 = [[],[(0,S1,1.0)],[(0,S1,2.0)],[]]
                # IntensiveMaximum => [[],[(0,S1,1.0/1.0)],[(0,S1,2.0/2.0)],[]]
                self.assertTrue(field.getArray().isEqual(DataArrayDouble([-1000.0, 1.0, 1.0, -1000.0]),1e-8))
                self.assertTrue( dec.retrieveNonFetchedIds().isEqual(DataArrayInt([0,3])) )
            # target -> source
            dec.sendData()

        # Some clean up that still needs MPI communication, so to be done before MPI_Finalize()
        dec.release()
        target_group.release()
        source_group.release()
        MPI.COMM_WORLD.Barrier()

def Source_Proc_0():
    coo = DataArrayDouble([(0,2),(2,2),(4,2),(0,4),(2,4),(4,4),(0,6),(2,6)])
    m = MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
    m.insertNextCell(NORM_QUAD4,[0,1,4,3])
    m.insertNextCell(NORM_QUAD4,[1,2,5,4])
    m.insertNextCell(NORM_QUAD4,[3,4,7,6])
    return m

def Source_Proc_1():
    coo = DataArrayDouble([(6,2),(8,2),(10,2),(6,4),(8,4),(10,4)])
    m = MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
    m.insertNextCell(NORM_QUAD4,[0,1,4,3])
    m.insertNextCell(NORM_QUAD4,[1,2,5,4])
    return m

def Target_Proc_2():
    coo = DataArrayDouble([(1,0),(3.5,0),(1,3),(3.5,3)])
    m = MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
    m.insertNextCell(NORM_QUAD4,[0,1,3,2])
    return m

def Target_Proc_3():
    coo = DataArrayDouble([(6,0),(7,0),(8,0),(9,0),(10,0),
                           (6,1),(7,1),(9,1),(10,1),
                           (7,3),(8,3),
                           (6,4),(7,4)])
    m = MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
    m.insertNextCell(NORM_QUAD4,[0,1,6,5])
    m.insertNextCell(NORM_QUAD4,[1,2,10,9])
    m.insertNextCell(NORM_QUAD4,[5,6,12,11])
    m.insertNextCell(NORM_QUAD4,[3,4,8,7])
    return m

if __name__ == "__main__":
    unittest.main()
    MPI.Finalize()

