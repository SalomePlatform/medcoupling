
.. _python_testMEDCouplingremapper1_solution:

Interpoler avec MEDCouplingRemapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import MEDCoupling as mc
	from MEDCouplingRemapper import MEDCouplingRemapper 
	# Target mesh
	arr = mc.DataArrayDouble(11)
	arr.iota(0)
	trgMesh = mc.MEDCouplingCMesh()
	trgMesh.setCoords(arr,arr)
	trgMesh = trgMesh.buildUnstructured()	
	# Source mesh
	arr = mc.DataArrayDouble(21)
	arr.iota(0)
	arr *= 0.5
	srcMesh = mc.MEDCouplingCMesh()
	srcMesh.setCoords(arr,arr)
	srcMesh = srcMesh.buildUnstructured()	
	# Triangularize some cells in source
	tmp = srcMesh[:20]    # Extract a sub-part of srcMesh
	tmp.simplexize(0)
	srcMesh = mc.MEDCouplingUMesh.MergeUMeshes([tmp,srcMesh[20:]])
	# Interpolate cells to cells
	remap = MEDCouplingRemapper()
	remap.prepare(srcMesh,trgMesh,"P0P0")
	# Check matrix
	myMatrix = remap.getCrudeMatrix()
	print myMatrix
	sumByRows = mc.DataArrayDouble(len(myMatrix))
	for i,wIt in enumerate(sumByRows):
	  su = 0.
	  for it in myMatrix[i]:
	    su += myMatrix[i][it]
	  wIt[0] = su
	print "Is interpolation well prepared?", sumByRows.isUniform(1.,1e-12)
	# Source field construction
	srcField = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
	srcField.setMesh(srcMesh)
	srcField.fillFromAnalytic(1,"7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))")
	srcField.getArray().setInfoOnComponent(0, "powercell [W]")
	# Transfer field
	#remap.transferField(srcField, 1e300)
	srcField.setNature(mc.IntensiveMaximum)
	trgFieldCV = remap.transferField(srcField,1e300)
	# IntensiveMaximum
	integSource = srcField.integral(True)[0]
	integTarget =  trgFieldCV.integral(True)[0]
	print "IntensiveMaximum -- integrals: %lf == %lf" % (integSource, integTarget)
	
	accSource = srcField.getArray().accumulate()[0]
	accTarget = trgFieldCV.getArray().accumulate()[0]
	print "IntensiveMaximum -- sums: %lf != %lf" % (accSource, accTarget)
	# ExtensiveConservation
	srcField.setNature(mc.ExtensiveConservation)
	trgFieldI = remap.transferField(srcField,1e300)
	#
	integSource = srcField.integral(True)[0]
	integTarget =  trgFieldI.integral(True)[0]
	print "ExtensiveConservation -- integrals: %lf != %lf" % (integSource, integTarget)
	
	accSource = srcField.getArray().accumulate()[0]
	accTarget = trgFieldI.getArray().accumulate()[0]
	print "ExtensiveConservation -- sums: %lf == %lf" % (accSource, accTarget)
