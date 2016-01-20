
.. _python_testMEDLoaderAdvancedAPI1_solution:

Reading, Writing a MED file using MEDLoader advanced API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import MEDLoader as ml
	# Mesh creation
	targetCoords = [-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
	targetConn = [0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
	targetMesh = ml.MEDCouplingUMesh("MyMesh",2)
	targetMesh.allocateCells(5)
	targetMesh.insertNextCell(ml.NORM_TRI3,3,targetConn[4:7])
	targetMesh.insertNextCell(ml.NORM_TRI3,3,targetConn[7:10])
	targetMesh.insertNextCell(ml.NORM_QUAD4,4,targetConn[0:4])
	targetMesh.insertNextCell(ml.NORM_QUAD4,4,targetConn[10:14])
	targetMesh.insertNextCell(ml.NORM_QUAD4,4,targetConn[14:18])
	myCoords = ml.DataArrayDouble(targetCoords,9,2)
	myCoords.setInfoOnComponents(["X [km]","YY [mm]"])
	targetMesh.setCoords(myCoords)
	# Build the 2D faces from the 3D volumes (descending connectivity)
	targetMeshConsti, _, _, _, _ = targetMesh.buildDescendingConnectivity()
	targetMesh1 = targetMeshConsti[[3,4,7,8]]
	targetMesh1.setName(targetMesh.getName())
	#
	# Meshes
	#
	meshMEDFile = ml.MEDFileUMesh()
	meshMEDFile.setMeshAtLevel(0,targetMesh)
	meshMEDFile.setMeshAtLevel(-1,targetMesh1)
	# Some groups on cells Level 0
	grp0_0 = ml.DataArrayInt([0,1,3]) 
	grp0_0.setName("grp0_Lev0")
	grp1_0 = ml.DataArrayInt([1,2,3,4])
	grp1_0.setName("grp1_Lev0")
	meshMEDFile.setGroupsAtLevel(0, [grp0_0,grp1_0])
	# Some groups on cells Level -1
	grp0_M1 = ml.DataArrayInt([0,1])
	grp0_M1.setName("grp0_LevM1")
	grp1_M1 = ml.DataArrayInt([0,1,2])
	grp1_M1.setName("grp1_LevM1")
	grp2_M1 = ml.DataArrayInt([1,2,3])
	grp2_M1.setName("grp2_LevM1")
	meshMEDFile.setGroupsAtLevel(-1,[grp0_M1,grp1_M1,grp2_M1])	
	# Write everything
	meshMEDFile.write("TargetMesh2.med",2) # 2 stands for write from scratch 
	# Re-read and test equality
	meshMEDFileRead = ml.MEDFileMesh.New("TargetMesh2.med")  # a new is needed because it returns a MEDFileUMesh (MEDFileMesh is abstract)
	meshRead0 = meshMEDFileRead.getMeshAtLevel(0)
	meshRead1 = meshMEDFileRead.getMeshAtLevel(-1)
	print "Is level 0 in the file equal to 'targetMesh'?", meshRead0.isEqual(targetMesh,1e-12)
	print "Is level 0 in the file equal to 'targetMesh1'?", meshRead1.isEqual(targetMesh1,1e-12)
	# Read groups
	print meshMEDFileRead.getGrpNonEmptyLevels("grp0_Lev0")
	grp0_0_read = meshMEDFileRead.getGroupArr(0,"grp0_Lev0")
	print "Is group 'grp0_Lev0' equal to what is read in the file?" , grp0_0_read.isEqual(grp0_0)
	#
	# Fields
	#
	f = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)
	f.setTime(5.6,7,8)
	f.setArray(targetMesh.computeCellCenterOfMass())
	f.setMesh(targetMesh)
	f.setName("AFieldName")
	# Prepare field for writing
	fMEDFile = ml.MEDFileField1TS()
	fMEDFile.setFieldNoProfileSBT(f)     # No profile desired on the field, Sort By Type
	# *Append* the field to an existing file
	fMEDFile.write("TargetMesh2.med",0) # 0 is very important here because we want to append to TargetMesh2.med and not to scratch it
	# Read the field
	fMEDFileRead = ml.MEDFileField1TS("TargetMesh2.med",f.getName(),7,8)
	fRead1 = fMEDFileRead.getFieldOnMeshAtLevel(ml.ON_CELLS,0,meshMEDFileRead) # Quickest way, not re-reading mesh in the file.
	fRead2 = fMEDFileRead.getFieldAtLevel(ml.ON_CELLS,0)                       # Like above, but this time the mesh is read!
	print "Does the field remain OK with the quick method?", fRead1.isEqual(f,1e-12,1e-12)
	print "Does the field remain OK with the slow method?", fRead2.isEqual(f,1e-12,1e-12)
	#
	# Writing and Reading fields on profile using MEDLoader advanced API
	#
	pfl = ml.DataArrayInt([1,2,3]) 
	pfl.setName("My1stPfl")
	fPart = f.buildSubPart(pfl)
	fPart.setName("fPart")
	#
	fMEDFile2 = ml.MEDFileField1TS()
	fMEDFile2.setFieldProfile(fPart,meshMEDFileRead,0,pfl) # 0 is the relative level (here 0 means 3D)
	fMEDFile2.write("TargetMesh2.med",0) # 0 is paramount to indicate that we *append* (and no overwrite) to the MED file
	#
	fMEDFileRead2 = ml.MEDFileField1TS("TargetMesh2.med",fPart.getName(),7,8)
	fPartRead, pflRead = fMEDFileRead2.getFieldWithProfile(ml.ON_CELLS,0,meshMEDFileRead)
	print "Is the partial field correclty read?", fPartRead.isEqualWithoutConsideringStr(fPart.getArray(),1e-12)
	print "Is the list of cell identifiers matching?", pflRead.isEqualWithoutConsideringStr(pfl)
	
