# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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

from MEDPartitioner import *
from MEDLoader import *
import unittest
from MEDLoaderDataForTest import MEDLoaderDataForTest

class MEDPartitionerTest(unittest.TestCase):
    def testPartition(self):
        fname="PyPartitionTest.med"
        data=MEDLoaderDataForTest.buildACompleteMEDDataStructureWithFieldsOnCells_1()
        data.write(fname,2)
        part_file=MEDPartitioner(fname,2)
        part_data=MEDPartitioner(data,2)
        part_file.write("splitted_PyPartitionTest_1")
        part_data.write("splitted_PyPartitionTest_2")
        part_file_xml=MEDPartitioner("splitted_PyPartitionTest_1.xml")
        part_data_xml=MEDPartitioner("splitted_PyPartitionTest_2.xml")
        data1=part_file_xml.getMEDFileData()
        data2=part_data_xml.getMEDFileData()
        m1d=data1.getMeshes().getMeshAtPos(0)
        m2d=data2.getMeshes().getMeshAtPos(0)
        self.assertTrue(m1d.isEqual(m2d,1e-12))
    pass
    def testPartitionGraph(self):
        data=MEDLoaderDataForTest.buildACompleteMEDDataStructureWithFieldsOnCells_1()
        m=data.getMeshes().getMeshAtPos(0)
        graph=MEDPartitioner.Graph(m.getLevel0Mesh().generateGraph())
        graph.partGraph(2)
        tool=MEDPartitioner(data,graph)
        data2=tool.getMEDFileData()
        self.assertEqual( 2, data2.getMeshes().getNumberOfMeshes() )
    pass
    def testPartitionWithJoints(self):
        # cartesian mesh 4x4
        arr=DataArrayDouble(5) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr)
        m=c.buildUnstructured()
        m.setName("mesh")
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        data=MEDFileData()
        data.setMeshes(ms)
        part_file=MEDPartitioner(data,4,"metis",True,True,True)
        data_file=part_file.getMEDFileData()
        meshes=data_file.getMeshes()
        self.assertEqual( 4, meshes.getNumberOfMeshes())
        self.assertEqual( 3, meshes.getMeshAtPos(0).getJoints().getNumberOfJoints())
        self.assertEqual( 3, meshes.getMeshAtPos(1).getJoints().getNumberOfJoints())
        self.assertEqual( 3, meshes.getMeshAtPos(2).getJoints().getNumberOfJoints())
        self.assertEqual( 3, meshes.getMeshAtPos(3).getJoints().getNumberOfJoints())
        joints=meshes.getMeshAtPos(0).getJoints()
        self.assertEqual( 1, joints.getJointAtPos(0).getDomainNumber(), 1)
        #VSR (10/05/2016): changed to work with metis 5.1... to be confirmed!
        #self.assertEqual( 2, joints.getJointAtPos(1).getDomainNumber(), 2)
        #self.assertEqual( 3, joints.getJointAtPos(2).getDomainNumber(), 3)
        self.assertEqual( 3, joints.getJointAtPos(1).getDomainNumber(), 3)
        self.assertEqual( 2, joints.getJointAtPos(2).getDomainNumber(), 2)
        self.assertEqual( 2, joints.getJointAtPos(0).getStepAtPos(0).getNumberOfCorrespondences())
        self.assertEqual( 2, joints.getJointAtPos(1).getStepAtPos(0).getNumberOfCorrespondences())
        self.assertEqual( 1, joints.getJointAtPos(2).getStepAtPos(0).getNumberOfCorrespondences())
        found=0
        for ii in range(joints.getJointAtPos(0).getStepAtPos(0).getNumberOfCorrespondences()):
            correspond=joints.getJointAtPos(0).getStepAtPos(0).getCorrespondenceAtPos(ii)
            #VSR (10/05/2016): changed to work with metis 5.1... to be confirmed!
            #if correspond.getCorrespondence().isEqual(DataArrayInt([1,3,2,4])):
            if correspond.getCorrespondence().isEqual(DataArrayInt([3,1,4,2])):
                found+=1
                self.assertEqual(NORM_QUAD4, correspond.getLocalGeometryType())
                self.assertEqual(NORM_QUAD4, correspond.getRemoteGeometryType())
            pass
        pass
        self.assertEqual(1,found)
    pass
    def testPartitionPartGraph(self):
        arr=DataArrayDouble(5) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr)
        m=c.buildUnstructured()
        part=MEDPartitioner.Graph(m.generateGraph())
        part.partGraph(2)
        a=part.getGraph()
        p=part.getPartition()
        self.assertTrue(isinstance(a,MEDCouplingSkyLineArray))
        self.assertTrue(isinstance(p,MEDCouplingSkyLineArray))
        self.assertTrue(part.nbVertices() > 0 )
    pass

if __name__ == "__main__":
  unittest.main()

