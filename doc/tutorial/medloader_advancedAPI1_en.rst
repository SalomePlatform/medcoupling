
Reading, Writing a MED file using MEDLoader's advanced API
----------------------------------------------------------

The advanced API is incarnated by the MEDFile* classes in the MEDLoader library.

* MEDFileMesh, MEDFileUMesh, MEDFileCMesh
* MEDFileMeshes, MEDFileMeshMultiTS
* MEDFileField1TS, MEDFileFieldMultiTS
* MEDFileFields, MEDFileFieldGlobs
* MEDFileData

Objective
~~~~~~~~~

Write a mesh and a field from scratch, re-read them and compare the result.

Topics covered:

* Read/Write Mesh using MEDLoader's advanced API
* Read/Write Field using MEDLoader's advanced API

Implementation start
~~~~~~~~~~~~~~~~~~~~

To implement this exercise we use the Python scripting language and import the MEDLoader Python module.
The whole MEDCoupling module is fully included in MEDLoader. No need to import MEDCoupling when MEDLoader has been loaded. ::

	import MEDLoader as ml

Writing and Reading meshes using MEDLoader's advanced API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First of all, creation of a mesh "targetMesh". ::

	targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        targetMesh=ml.MEDCouplingUMesh.New("MyMesh",2)
        targetMesh.allocateCells(5)
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
	targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        myCoords=ml.DataArrayDouble.New(targetCoords,9,2)
        targetMesh.setCoords(myCoords)
        

.. note:: targetMesh is grouped by geometric type.

Build "targetMesh1" representing the sub-constituents (faces) of "targetMesh" reduced to cell ids [3,4,7,8]. 
::

	targetMeshConsti=targetMesh.buildDescendingConnectivity()[0]
	targetMesh1=targetMeshConsti[[3,4,7,8]]
	targetMesh1.setName(targetMesh.getName())

.. note:: "targetMesh1" will be recorded as a part of the same global mesh in the MED file, so it must have the same name!

Then we are ready to write targetMesh and targetMesh1 into TargetMesh2.med. ::

	meshMEDFile=MEDFileUMesh.New()
	meshMEDFile.setMeshAtLevel(0,targetMesh)
	meshMEDFile.setMeshAtLevel(-1,targetMesh1)
	meshMEDFile.write("TargetMesh2.med",2) # 2 stands for write from scratch

Create 2 groups on level 0. The first called "grp0_Lev0" on cells [0,1,3] and the second called "grp1_Lev0" on cells [1,2,3,4] ::	

	grp0_0=ml.DataArrayInt.New([0,1,3]) ; grp0_0.setName("grp0_Lev0")
	grp1_0=ml.DataArrayInt.New([1,2,3,4]) ; grp1_0.setName("grp1_Lev0")
	meshMEDFile.setGroupsAtLevel(0,[grp0_0,grp1_0])

Create 3 groups on level -1. The 1st called "grp0_LevM1" on cells [0,1], the 2nd called "grp1_LevM1" on cells [0,1,2], and the 3rd called "grp2_LevM1" on cells [1,2,3] ::

	grp0_M1=ml.DataArrayInt.New([0,1]) ; grp0_M1.setName("grp0_LevM1")
	grp1_M1=ml.DataArrayInt.New([0,1,2]) ; grp1_M1.setName("grp1_LevM1")
	grp2_M1=ml.DataArrayInt.New([1,2,3]) ; grp2_M1.setName("grp2_LevM1")
	meshMEDFile.setGroupsAtLevel(-1,[grp0_M1,grp1_M1,grp2_M1])
	

Then trying to read it. ::

	meshMEDFileRead=MEDFileMesh.New("TargetMesh2.med")
	meshRead0=meshMEDFileRead.getMeshAtLevel(0)
	meshRead1=meshMEDFileRead.getMeshAtLevel(-1)
	print "Is the mesh at level 0 read in file equals targetMesh ? %s"%(meshRead0.isEqual(targetMesh,1e-12))
	print "Is the mesh at level -1 read in file equals targetMesh ? %s"%(meshRead1.isEqual(targetMesh1,1e-12))

Print available levels for group "grp0_Lev0" ::

	print meshMEDFileRead.getGrpNonEmptyLevels("grp0_Lev0")

Request for cell ids of group "grp0_Lev0" ::

	grp0_0_read=meshMEDFileRead.getGroupArr(0,"grp0_Lev0")
	print "Is group \"grp0_Lev0\" are the same ? %s"%(grp0_0_read.isEqual(grp0_0))

Writing and Reading fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

Creation of a simple vector field on cells called f.  ::

	f=ml.MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
	f.setTime(5.6,7,8)
	f.setArray(targetMesh.computeCellCenterOfMass())
	f.setMesh(targetMesh)
	f.setName("AFieldName")

Put f into a MEDFileField1TS for preparation of MED writing ::

	fMEDFile=MEDFileField1TS.New()
	fMEDFile.setFieldNoProfileSBT(f)

Append field to "TargetMesh2.med" ::

	fMEDFile.write("TargetMesh2.med",0) # 0 is very important here because we want to append to TargetMesh2.med and not to overwrite it

Read it : ::

	fMEDFileRead=MEDFileField1TS.New("TargetMesh2.med",f.getName(),7,8)
	fRead1=fMEDFileRead.getFieldOnMeshAtLevel(ON_CELLS,0,meshMEDFileRead) # fastest method. No reading of the supporting mesh.
	fRead2=fMEDFileRead.getFieldAtLevel(ON_CELLS,0) # like above but mesh is re-read from file...
	print "Does the field f remain the same using fast method ? %s"%(fRead1.isEqual(f,1e-12,1e-12))
	print "Does the field f remain the same using slow method ? %s"%(fRead2.isEqual(f,1e-12,1e-12))
	
Writing and Reading fields on a "profile"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Build a reduction on cells [1,2,3] of f and call it fPart. ::

	pfl=ml.DataArrayInt.New([1,2,3]) ; pfl.setName("My1stPfl")
	fPart=f.buildSubPart(pfl)
	fPart.setName("fPart")

Put it into MEDFileField1TS data structure. ::

	fMEDFile2=MEDFileField1TS.New()
	fMEDFile2.setFieldProfile(fPart,meshMEDFileRead,0,pfl)
	fMEDFile2.write("TargetMesh2.med",0) # 0 is very important here because we want to append to TargetMesh2.med and not to scratch it

Read "fPart" field from File "TargetMesh2.med". ::

	fMEDFileRead2=MEDFileField1TS.New("TargetMesh2.med",fPart.getName(),7,8)
	fPartRead,pflRead=fMEDFileRead2.getFieldWithProfile(ON_CELLS,0,meshMEDFileRead)
	print fPartRead.isEqualWithoutConsideringStr(fPart.getArray(),1e-12)
	print pflRead.isEqualWithoutConsideringStr(pfl)

Solution
~~~~~~~~

:ref:`python_testMEDLoaderAdvancedAPI1_solution`
