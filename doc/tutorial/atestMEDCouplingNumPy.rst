
.. _python_testMEDCouplingNumPy_solution:

Playing with NumPy and SciPy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
	
	import MEDCoupling as mc
	
	#
	# NumPy
	#
	import numpy as np
	
	# Checking NumPy binding
	assert(mc.MEDCouplingHasNumPyBindings())
	# Playing with conversion and shared data
	arr = mc.DataArrayDouble(12)
	arr[:] = 4.
	nparr = arr.toNumPyArray()
	print nparr.__repr__()
	print nparr.tolist()
	nparr[::2] = 7.
	print nparr.__repr__()
	print arr.__repr__()
	del arr
	import gc; gc.collect()     # Make sure the object has been deleted
	print nparr.__repr__()
	arr2 = mc.DataArrayDouble(nparr)
	print arr2.__repr__()
	nparr[:] = 5.
	print nparr.__repr__()
	print arr2.__repr__()
	# Writing to file
	f = open("toto.data","w+b")
	a = np.memmap(f,dtype='float64',mode='w+',offset=0,shape=nparr.shape)
	a[:] = nparr[:]
	f.flush()
	# Re-reading file
	f2 = open("toto.data","r+b")
	b = np.memmap(f2,dtype='float64',mode='r',offset=0,shape=(12,))
	a[:] = 3.14
	f.flush()
	b = np.memmap(f2,dtype='float64',mode='r',offset=0,shape=(12,))
	print b.__repr__()
	#
	# SciPy
	#
	assert(mc.MEDCouplingHasSciPyBindings())
	c1 = mc.MEDCouplingCMesh()
	arr1 = mc.DataArrayDouble(7) 
	arr1.iota() 
	c1.setCoords(arr1,arr1,arr1)
	c2 = mc.MEDCouplingCMesh()
	arr2 = mc.DataArrayDouble(9)
	arr2.iota() 
	arr2 *= 6./8.
	c2.setCoords(arr2,arr2,arr2)
	c1 = c1.buildUnstructured()
	c2 = c2.buildUnstructured()
	c2.translate([6.,0.,0.])
	c = mc.MEDCouplingUMesh.MergeUMeshes([c1,c2])
	c.mergeNodes(1e-12)
	skinAndNCFaces = c.computeSkin()
	skinAndNCFaces.zipCoords()
	# Isolating non conform cells
	from MEDCouplingRemapper import MEDCouplingRemapper
	rem = MEDCouplingRemapper()
	rem.setMaxDistance3DSurfIntersect(1e-12)
	rem.setMinDotBtwPlane3DSurfIntersect(0.99)
	rem.prepare(skinAndNCFaces,skinAndNCFaces,"P0P0")
	mat = rem.getCrudeCSRMatrix()
	indptr = mc.DataArrayInt(mat.indptr)
	indptr2 = indptr.deltaShiftIndex()
	cellIdsOfSkin = indptr2.findIdsEqual(1)
	skin = skinAndNCFaces[cellIdsOfSkin]
	skin.writeVTK("skin.vtu")
