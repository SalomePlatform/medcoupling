
.. _python_testMEDCouplingfielddouble1_solution:

Playing with fields
~~~~~~~~~~~~~~~~~~~

::

	import MEDCoupling as mc
	
	# Create an unstructured mesh from a Cartesian one
	xarr = mc.DataArrayDouble.New(11,1)
	xarr.iota(0.)
	cmesh = mc.MEDCouplingCMesh.New()
	cmesh.setCoords(xarr,xarr,xarr)
	mesh = cmesh.buildUnstructured()
	mesh.convertToPolyTypes(mc.DataArrayInt.Range(0,mesh.getNumberOfCells(),2))
	# Create a field
	f = mesh.fillFromAnalytic(mc.ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")  # 1 means that the field should have one component
	f.setName("MyField")
	# A variant: 
	f2 = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
	f2.setMesh(mesh)
	f2.setName("MyField2")
	f2.fillFromAnalytic(1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")    # 1 means that the field should have one component
	print "Are f and f2 equal?", f.isEqualWithoutConsideringStr(f2,1e-12,1e-12)
	#
	da1 = f.getArray()              # a DataArrayDouble, which is a direct reference (not a copy) of the field's values
	ids1 = da1.findIdsInRange(0.,5.)
	fPart1 = f.buildSubPart(ids1)
	fPart1.writeVTK("ExoField_fPart1.vtu")
	ids2 = f.getArray().findIdsInRange(50.,1.e300)
	fPart2 = f.buildSubPart(ids2)
	# Renumbering cells to follow MED file rules
	fPart1Cpy = fPart1.deepCopy()
	o2n = fPart1Cpy.getMesh().sortCellsInMEDFileFrmt()
	fPart1Cpy.getArray().renumberInPlace(o2n)
	# Check that fPart1Cpy and fPart1 are the same
	fPart1Cpy.substractInPlaceDM(fPart1,12,1e-12)
	fPart1Cpy.getArray().abs()
	print "Are the fields equal?", (fPart1Cpy.getArray().accumulate()[0]<1e-12)
	# Aggregate fields
	fPart12 = mc.MEDCouplingFieldDouble.MergeFields([fPart1,fPart2])
	fPart12.writeVTK("ExoField_fPart12.vtu")
	# Evaluation on points
	bary = fPart12.getMesh().computeCellCenterOfMass()
	arr1 = fPart12.getValueOnMulti(bary)
	arr2 = f.getValueOnMulti(bary)
	delta = arr1-arr2
	delta.abs()
	print "Is field evaluation matching?", (delta.accumulate()[0]<1e-12)
	# ExtensiveMaximum computations
	integ1 = fPart12.integral(0,True)
	integ1_bis = fPart12.getArray().accumulate()[0]
	print "First integral matching ?", ( abs(integ1 - integ1_bis) < 1e-8 )
	fPart12.getMesh().scale([0.,0.,0.], 1.2)	
	integ2 = fPart12.integral(0,True)
	print "Second integral matching ?", ( abs(integ2-integ1_bis*1.2*1.2*1.2) < 1e-8 )
	# Explosion of field
	fVec = mesh.fillFromAnalytic(mc.ON_CELLS,3,"(x-5.)*IVec+(y-5.)*JVec+(z-5.)*KVec")
	fVecPart1 = fVec.buildSubPart(ids1)
	fVecPart1.setName("fVecPart1")
	cells = fPart1.getMesh().getNumberOfCells() * [None]
	for icell,vec in enumerate(fVecPart1.getArray()):
	  m = fPart1.getMesh()[[icell]]
	  m.zipCoords()      # Not mandatory but saves memory
	  m.translate(vec)
	  cells[icell] = m
	  pass
	meshFVecPart1Exploded = mc.MEDCouplingUMesh.MergeUMeshes(cells)
	fPart1.setMesh(meshFVecPart1Exploded)
	fPart1.writeVTK("ExoField_fPart1_explo.vtu")
	

