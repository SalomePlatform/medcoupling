
.. _python_testMEDLoaderBasicAPI1_solution:

Reading, Writing a MED file using MEDLoader basic API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import medcoupling as mc
	# Mesh creation
	targetCoords = [-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
	targetConn = [0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
	targetMesh = mc.MEDCouplingUMesh("MyMesh",2)
	targetMesh.allocateCells(5)
	targetMesh.insertNextCell(mc.NORM_TRI3,3,targetConn[4:7])
	targetMesh.insertNextCell(mc.NORM_TRI3,3,targetConn[7:10])
	targetMesh.insertNextCell(mc.NORM_QUAD4,4,targetConn[0:4])
	targetMesh.insertNextCell(mc.NORM_QUAD4,4,targetConn[10:14])
	targetMesh.insertNextCell(mc.NORM_QUAD4,4,targetConn[14:18])
	myCoords = mc.DataArrayDouble(targetCoords,9,2)
	myCoords.setInfoOnComponents(["X [km]","YY [mm]"])
	targetMesh.setCoords(myCoords)
	# Writing mesh only
	mc.WriteUMesh("TargetMesh.med",targetMesh,True)  # True means 'from scratch'
	# Re-read it and test equality
	meshRead = mc.ReadUMeshFromFile("TargetMesh.med",targetMesh.getName(),0)
	print("Is the read mesh equal to 'targetMesh' ?", meshRead.isEqual(targetMesh,1e-12))
	# Writing a field and its support mesh in one go
	f = mc.MEDCouplingFieldDouble.New(mc.ON_CELLS, mc.ONE_TIME)
	f.setTime(5.6,7,8)                              # Declare the timestep associated to the field
	f.setArray(targetMesh.computeCellCenterOfMass())
	f.setMesh(targetMesh)
	f.setName("AFieldName")
	mc.WriteField("MyFirstField.med",f,True)
	# Re-read it and test equality
	f2 = mc.ReadFieldCell("MyFirstField.med", f.getMesh().getName(), 0, f.getName(), 7, 8)
	print("Is the read field identical to 'f' ?", f2.isEqual(f,1e-12,1e-12))
	# Writing in several steps
	mc.WriteUMesh("MySecondField.med",f.getMesh(),True)
	mc.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f)
	# A second field to write
	f2 = f.clone(True)         # 'True' means that we need a deep copy
	f2.getArray()[:] = 2.0
	f2.setTime(7.8,9,10)
	mc.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f2)
	# Re-read and test this two-timestep field
	f3 = mc.ReadFieldCell("MySecondField.med",f.getMesh().getName(),0,f.getName(),7,8)
	print("Is the field read in file equals to 'f' ?", f.isEqual(f3,1e-12,1e-12))
	f4 = mc.ReadFieldCell("MySecondField.med",f.getMesh().getName(),0,f.getName(),9,10)
	print("Is the field read in file equals to 'f2' ?", f2.isEqual(f4,1e-12,1e-12))

