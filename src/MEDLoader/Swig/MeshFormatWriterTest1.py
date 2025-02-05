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

from MEDLoader import *
import unittest
from MEDLoaderDataForTest import MEDLoaderDataForTest,WriteInTmpDir

class MeshFormatWriterTest2(unittest.TestCase):

  def createMesh(self, nb_seg):

    # Create mesh
    mesh_dim = 3
    self.bounding_box_coords = [(0, 100)]*mesh_dim
    # with 25 segs per side, time of test run
    # 9.12: 7.8 s
    # 9.13: 0.1 s

    # with 100 nb_segs per side:
    # 9.12: more than 3 h
    # 9.13: 8.14 s

    nb_segs = [nb_seg]*mesh_dim
    box_sizes = [size_coords[1]-size_coords[0] for size_coords in self.bounding_box_coords]
    print("box_sizes: ", box_sizes)
    box_steps = [box_sizes[i]/nb_segs[i] for i in range(mesh_dim)]
    l_nb_nodes = [nb_segs_i+1 for nb_segs_i in nb_segs]
    l_start = [size_coords[0] for size_coords in self.bounding_box_coords]
    l_end = [size_coords[1] for size_coords in self.bounding_box_coords]

    mesh = MEDCouplingIMesh("Mesh", mesh_dim, l_nb_nodes, l_start, box_steps)
    self.mesh = mesh.buildUnstructured()
    pass

  def createVolumeGroups(self):
    nb_cells = self.mesh.getNumberOfCells()
    nb_half_cells = int(round(nb_cells/2))
    grp0 = DataArrayInt(list(range(nb_half_cells)))
    grp0.setName("half1")
    grp1 = DataArrayInt(list(range(nb_half_cells, nb_cells)))
    grp1.setName("half2")
    self.volumeGroups = [grp0, grp1]
    pass

  def createSkinGroups(self):
    # create groups on each side of skin mesh
    nb_faces = self.skinMesh.getNumberOfCells()
    #xMin, xMax, yMin, yMax, zMin, zMax
    xMin, xMax = self.bounding_box_coords[0]
    yMin, yMax = self.bounding_box_coords[1]
    zMin, zMax = self.bounding_box_coords[2]
    tol = 1e-10
    barycenters = self.skinMesh.computeIsoBarycenterOfNodesPerCell()
    left = []
    right = []
    bottom = []
    top = []
    back = []
    front = []
    for i in range(nb_faces):
      coord = barycenters[i]
      x, y, z = coord.getValues()
      if abs(x-xMin) < tol:
        left.append(i)
      elif abs(x-xMax) < tol:
        right.append(i)
      elif abs(y-yMin) < tol:
        back.append(i)
      elif abs(y-yMax) < tol:
        front.append(i)
      elif abs(z-zMin) < tol:
        bottom.append(i)
      elif abs(z-zMax) < tol:
        top.append(i)

    grp_left = DataArrayInt(left)
    grp_left.setName("left")
    grp_right = DataArrayInt(right)
    grp_right.setName("right")
    grp_back = DataArrayInt(back)
    grp_back.setName("back")
    grp_front = DataArrayInt(front)
    grp_front.setName("front")
    grp_bottom = DataArrayInt(bottom)
    grp_bottom.setName("bottom")
    grp_top = DataArrayInt(top)
    grp_top.setName("top")
    self.skinGroups = [grp_left, grp_right, grp_back, grp_front, grp_bottom, grp_top]
    pass

  def createField(self):
    # create a field on nodes with 1 component
    coords = self.mesh.getCoords()
    coords_z = coords[:,2]
    zMin, zMax = self.bounding_box_coords[2]
    coords_z_normed = coords_z/(zMax-zMin)
    size_min = 1e-4
    size_max = 1e-3
    coords_z_normed.applyFuncOnThis('(%f-%f)*tanh(x)/tanh(1)+%f'%(size_max, size_min, size_min))

    # Create the field with this array
    self.field = MEDCouplingFieldDouble(ON_NODES,ONE_TIME)
    self.field.setName("Elevation")
    self.field.setMesh(self.mesh)
    self.field.setTime(0,1,1)
    self.field.setArray(coords_z_normed)
    pass

  def createMEDMeshFile(self):
    self.meshMEDFile = MEDFileUMesh()
    self.meshMEDFile.setMeshAtLevel(0, self.mesh)
    self.meshMEDFile.setGroupsAtLevel(0, self.volumeGroups)
    self.meshMEDFile.setMeshAtLevel(-1, self.skinMesh)
    self.meshMEDFile.setGroupsAtLevel(-1, self.skinGroups)
    pass

  def setMeshInMEDFileData(self):
    # Set mesh in file data structure
    ms=MEDFileMeshes()
    ms.pushMesh(self.meshMEDFile)
    self.medFileData = MEDFileData()
    self.medFileData.setMeshes(ms)
    pass

  def setFieldInMEDFileData(self):
    mf=MEDFileFields()
    fmts0=MEDFileFieldMultiTS()
    mf.pushField(fmts0)
    f1ts=MEDFileField1TS()
    f1ts.setFieldNoProfileSBT(self.field)
    fmts0.pushBackTimeStep(f1ts)
    self.medFileData.setFields(mf)
    pass

  def writeToMeshFile(self, meshFileName):
    tmpMeshWriter = MeshFormatWriter()
    tmpMeshWriter.setMeshFileName(meshFileName)
    tmpMeshWriter.setMEDFileDS(self.medFileData)
    tmpMeshWriter.write()
    pass

  def readMeshFile(self, meshFileName):
    meshFormatReader = MeshFormatReader()
    meshFormatReader.setFile(meshFileName)
    medFileData = meshFormatReader.loadInMedFileDS()
    meshes = medFileData.getMeshes()
    meshNames = meshes.getMeshesNames()
    assert len(meshNames) == 1
    self.meshMEDFileOut = meshes.getMeshWithName(meshNames[0])

  def compareFamilies(self, levels):
    meshRef = self.meshMEDFile
    meshRead = self.meshMEDFileOut
    # compare families on levels
    for level in levels:
      print("level: ", level)
      # ids of families at this level
      referenceFamilyIds = meshRef.getFamilyFieldAtLevel(level).getDifferentValues().getValues()
      print("referenceFamilyIds; ", referenceFamilyIds)
      readFamilyIds      = meshRead.getFamilyFieldAtLevel(level).getDifferentValues().getValues()
      print("readFamilyIds: ", readFamilyIds)
      assert len(readFamilyIds) == len(readFamilyIds)

      for ref_id in referenceFamilyIds:
        # by convention, the id in the .mesh file is the abs of the family
        read_id = abs(ref_id)
        referenceName = meshRef.getFamilyNameGivenId(ref_id)
        readName = meshRead.getFamilyNameGivenId(read_id)

        referenceFamilyArr = meshRef.getFamilyArr(level, referenceName)
        readFamilyArr = meshRead.getFamilyArr(level, readName)

        # compare family arrays
        referenceFamilyArr = meshRef.getFamilyArr(level, referenceName)
        readFamilyArr = meshRead.getFamilyArr(level, readName)
        referenceFamilyArr.isEqualWithoutConsideringStr(readFamilyArr)
        #print(readFamilyArr_i)

        # compare family meshes
        referenceFamily = meshRef.getFamily(level, referenceName)
        readFamily = meshRead.getFamily(level, readName)

        assert referenceFamily.isEqualWithoutConsideringStr(readFamily, 1e-12)

  @WriteInTmpDir
  def testMeshHexaWithoutField(self):
    """
    Test writing .mesh without field does not crash (as it did before 9.13)
    """

    nb_seg = 3
    self.createMesh(nb_seg)

    # Create groups
    self.createVolumeGroups()
    # Build 1D mesh and create groups of faces
    self.skinMesh = self.mesh.computeSkin()
    self.createSkinGroups()
    self.createField()
    self.createMEDMeshFile()

    self.setMeshInMEDFileData()

    # write to med (to debug if needed)
    #medFileData.write("Mesh3D_hexa.med", 2)

    # write to .mesh
    self.writeToMeshFile("Mesh3D_hexa_%i.mesh"%nb_seg)
    # test succeeds if there is no crash
    pass

  @WriteInTmpDir
  def testMeshHexaWithField(self):
    """
    Test writing .mesh is not too slow (more than 9 seconds before 9.13)
    """

    nb_seg = 25
    self.createMesh(nb_seg)

    # Create groups
    self.createVolumeGroups()
    # Build 1D mesh and create groups of faces
    self.skinMesh = self.mesh.computeSkin()
    self.createSkinGroups()
    self.createField()
    self.createMEDMeshFile()

    # Set mesh in file data structure
    self.setMeshInMEDFileData()
    # Set field in file data structure
    self.setFieldInMEDFileData()

    # Write to med (to debug if needed)
    #medFileData.write("Mesh3D_hexa.med", 2)

    # Write to .mesh
    meshFileName = "Mesh3D_hexa_%i.mesh"%nb_seg
    self.writeToMeshFile(meshFileName)

    # Read Mesh File
    self.readMeshFile(meshFileName)

    # Compare families at levels
    self.compareFamilies([0, -1])

    pass

  @WriteInTmpDir
  def testMeshTetraWithField(self):
    """
    Test writing .mesh tetra
    """

    nb_seg = 3
    self.createMesh(nb_seg)

    self.mesh = self.mesh.tetrahedrize(PLANAR_FACE_5)[0].buildUnstructured()

    # Create groups
    self.createVolumeGroups()
    # Build 1D mesh and create groups of faces
    self.skinMesh = self.mesh.computeSkin()
    self.createSkinGroups()
    self.createField()
    self.createMEDMeshFile()

    # Set mesh in file data structure
    self.setMeshInMEDFileData()
    # Set field in file data structure
    self.setFieldInMEDFileData()

    # Write to med (to debug if needed)
    #medFileData.write("Mesh3D_tetra.med", 2)

    # Write to .mesh
    meshFileName = "Mesh3D_hexa_%i.mesh"%nb_seg
    self.writeToMeshFile(meshFileName)

    # Read Mesh File
    self.readMeshFile(meshFileName)

    # Compare families at levels
    self.compareFamilies([0, -1])

    pass

if __name__ == '__main__':
  unittest.main()
