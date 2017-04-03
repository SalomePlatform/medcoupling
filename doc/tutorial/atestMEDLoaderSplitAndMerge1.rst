
.. _python_testMEDLoaderSplitAndMerge1_solution:

Splitting and Merging a MED file using MEDLoader
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import MEDLoader as ml
	
	m0 = ml.MEDCouplingCMesh()
	arr = ml.DataArrayDouble(31,1) ; arr.iota(0.)
	m0.setCoords(arr,arr)
	m0 = m0.buildUnstructured()
	m00 = m0[::2]      # Extract even cells
	m00.simplexize(0) 
	m01 = m0[1::2]
	m0 = ml.MEDCouplingUMesh.MergeUMeshes([m00,m01])
	m0.getCoords()[:] *= 1/15.
	m0.setName("mesh")
	# Cell field
	cellField = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME) 
	cellField.setTime(5.6,5,6)
	cellField.setMesh(m0)
	cellField.setName("CellField")
	cellField.fillFromAnalytic(1,"exp(-((x-1)*(x-1)+(y-1)*(y-1)))")
	cellField.getArray().setInfoOnComponent(0,"powercell [W]")
	# Node field
	nodeField = ml.MEDCouplingFieldDouble(ml.ON_NODES,ml.ONE_TIME) 
	nodeField.setTime(5.6,5,6)
	nodeField.setMesh(m0)
	nodeField.setName("NodeField")
	nodeField.fillFromAnalytic(1,"exp(-((x-1)*(x-1)+(y-1)*(y-1)))")
	nodeField.getArray().setInfoOnComponent(0,"powernode [W]")
	# Splitting
	proc0 = m0.getCellsInBoundingBox([(0.,0.4),(0.,0.4)],1e-10)
	proc1 = proc0.buildComplement(m0.getNumberOfCells())
	#
	nodeField0 = nodeField[proc0] ; cellField0 = cellField[proc0] ; cellField0.setMesh(nodeField0.getMesh())
	nodeField1 = nodeField[proc1] ; cellField1 = cellField[proc1] ; cellField1.setMesh(nodeField1.getMesh())
	
	proc0_fname = "proc0.med"
	ml.WriteField(proc0_fname, nodeField0, True)
	ml.WriteFieldUsingAlreadyWrittenMesh(proc0_fname, cellField0)
	
	proc1_fname = "proc1.med"
	ml.WriteField(proc1_fname,nodeField1,True)
	ml.WriteFieldUsingAlreadyWrittenMesh(proc1_fname,cellField1)
	#
	# Merging - Sub-optimal method
	#
	cellField0_read = ml.ReadFieldCell("proc0.med","mesh",0,"CellField",5,6)
	cellField1_read = ml.ReadFieldCell("proc1.med","mesh",0,"CellField",5,6)
	cellField_read = ml.MEDCouplingFieldDouble.MergeFields([cellField0_read,cellField1_read])
	cellFieldCpy = cellField.deepCopy()
	cellFieldCpy.substractInPlaceDM(cellField_read,10,1e-12)
	cellFieldCpy.getArray().abs()
	print cellFieldCpy.getArray().isUniform(0.,1e-12)
	#
	nodeField0_read = ml.ReadFieldNode("proc0.med","mesh",0,"NodeField",5,6)
	nodeField1_read = ml.ReadFieldNode("proc1.med","mesh",0,"NodeField",5,6)
	nodeField_read = ml.MEDCouplingFieldDouble.MergeFields([nodeField0_read, nodeField1_read])
	nodeField_read.mergeNodes(1e-10)
	nodeFieldCpy = nodeField.deepCopy()
	nodeFieldCpy.mergeNodes(1e-10)
	nodeFieldCpy.substractInPlaceDM(nodeField_read,10,1e-12)
	print nodeFieldCpy.getArray().isUniform(0.,1e-12)
	#
	# Merging - Optimal method
	#
	fileNames = ["proc0.med","proc1.med"]
	msML = [ml.MEDFileMesh.New(fname) for fname in fileNames]
	fsML = [ml.MEDFileFields.New(fname) for fname in fileNames]
	mergeMLMesh = ml.MEDFileUMesh()
	mergeMLFields = ml.MEDFileFields()
	for lev in msML[0].getNonEmptyLevels():
		o2nML = len(msML[0].getNonEmptyLevels())*[None]
		cs = [mML.getCoords() for mML in msML]
		mergeMLMesh.setCoords(ml.DataArrayDouble.Aggregate(cs))
		ms = [mML.getMeshAtLevel(lev) for mML in msML]
		m = ml.MEDCouplingUMesh.MergeUMeshes(ms) ; m.setCoords(mergeMLMesh.getCoords())
		o2nML[lev] = m.sortCellsInMEDFileFrmt()
		mergeMLMesh.setMeshAtLevel(lev,m)
		pass
	
	for fieldName in fsML[0].getFieldsNames():
		fmts = [fML[fieldName] for fML in fsML]
		mergeField = ml.MEDFileFieldMultiTS()
		for dt,it,tim in fmts[0].getTimeSteps():
			fts = [fmt[dt,it] for fmt in fmts]
			arrs = len(fts)*[None]
			for typp in fts[0].getTypesOfFieldAvailable():
				arr1s = []
				if typp == ml.ON_CELLS:
					for ft in fts:
						for geoTyp,smth in ft.getFieldSplitedByType():
							if geoTyp != ml.NORM_ERROR:
								smth1 = filter(lambda x:x[0] == ml.ON_CELLS,smth)
								arr2s = [ft.getUndergroundDataArray()[elt[1][0]:elt[1][1]] for elt in smth1]
								arr1s.append(ml.DataArrayDouble.Aggregate(arr2s))
								pass
							pass
						pass
					pass
				else:
					for ft in fts:
						smth = filter(lambda x:x[0] == ml.NORM_ERROR,ft.getFieldSplitedByType())
						arr2 = ml.DataArrayDouble.Aggregate([ft.getUndergroundDataArray()[elt[1][0][1][0]:elt[1][0][1][1]] for elt in smth])
						arr1s.append(arr2)
						pass
					pass
				arr = ml.DataArrayDouble.Aggregate(arr1s)
				if typp == ml.ON_CELLS:
				     arr.renumberInPlace(o2nML[lev])
				mcf = ml.MEDCouplingFieldDouble(typp,ml.ONE_TIME) ; mcf.setName(fieldName) ; mcf.setTime(tim,dt,it) ; mcf.setArray(arr)
				mcf.setMesh(mergeMLMesh.getMeshAtLevel(lev)) ; mcf.checkConsistencyLight()
				mergeField.appendFieldNoProfileSBT(mcf)
				pass
			pass
		mergeMLFields.pushField(mergeField)
		pass
	mergeMLMesh.write("merge.med",2)
	mergeMLFields.write("merge.med",0)
