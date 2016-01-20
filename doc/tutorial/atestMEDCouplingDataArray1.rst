
.. _python_testMEDCouplingdataarray1_solution:

Playing with regular hexagons using DataArrayDouble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import MEDCoupling as mc 
	import math
	# Building the coordinates of the initial hexagon, centered at 0,0
	d = mc.DataArrayDouble(6,2)
	d[:,0] = 3.
	d[:,1] = range(6)
	d[:,1] *= math.pi/3.
	d = d.fromPolarToCart()
	d.setInfoOnComponents(["X [m]","Y [m]"])
	print d.getValues()
	print d
	print "Uniform array?", d.magnitude().isUniform(3.,1e-12)
	# Translating the 7 hexagons with a translation
	radius = 3.
	translationToPerform = [[0.,0.],[3./2.*radius,-radius*math.sqrt(3.)/2],[3./2.*radius,radius*math.sqrt(3.)/2],[0.,radius*math.sqrt(3.)],[-3./2.*radius,radius*math.sqrt(3.)/2],[-3./2.*radius,-radius*math.sqrt(3.)/2],[0.,-radius*math.sqrt(3.)]]
	ds = len(translationToPerform)*[None]
	for pos,t in enumerate(translationToPerform):
			 ds[pos] = d[:]		# Perform a deep copy of d and place it at position 'pos' in ds
			 ds[pos] += t		  # Adding a vector to a set of coordinates does a translation
			 pass
	# Identifying duplicate tuples
	d2 = mc.DataArrayDouble.Aggregate(ds)
	oldNbOfTuples = d2.getNumberOfTuples()
	c,cI = d2.findCommonTuples(1e-12)
	tmp = c[cI[0]:cI[0+1]]
	print tmp
	a = cI.deltaShiftIndex()
	b = a - 1
	myNewNbOfTuples = oldNbOfTuples - sum(b.getValues())
	o2n, newNbOfTuples = mc.DataArrayInt.ConvertIndexArrayToO2N(oldNbOfTuples,c,cI)
	print "Have I got the right number of tuples?"
	print "myNewNbOfTuples = %d, newNbOfTuples = %d" % (myNewNbOfTuples, newNbOfTuples)
	assert(myNewNbOfTuples == newNbOfTuples)
	# Extracting the unique set of tuples 
	d3 = d2.renumberAndReduce(o2n, newNbOfTuples)
	n2o = o2n.invertArrayO2N2N2O(newNbOfTuples)
	d3_bis = d2[n2o]
	print "Are d3 and d3_bis equal ? %s" % (str(d3.isEqual(d3_bis, 1e-12)))
	# Now translate everything
	d3 += [3.3,4.4]
	# And build an unstructured mesh representing the final pattern
	m = mc.MEDCouplingUMesh("My7hexagons",2)
	m.setCoords(d3)
	print "Mesh dimension is", m.getMeshDimension()
	print "Spatial dimension is", m.getCoords().getNumberOfComponents()
	m.allocateCells(7)
	for i in xrange(7):
		cell_connec = o2n[6*i:6*(i+1)]
		m.insertNextCell(mc.NORM_POLYGON, cell_connec.getValues())
		pass
	# Check that everything is coherent (will throw if not)
	m.checkConsistencyLight()
	# Write the result into a VTU file that can be read with ParaView
	m.writeVTK("My7hexagons.vtu")

