
.. _python_testmedcouplingloaderex2_solution:

Intersection géométrique de maillages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import MEDLoader as ml
	
	def displayVTK(m,fname):
		tmp = m.deepCopy()
		tmp.tessellate2D(0.1)
		tmp.writeVTK(fname)
		return

	# Read and clean Fixe.med
	fixe = ml.MEDFileMesh.New("Fixe.med")
	fixm = fixe.getMeshAtLevel(0)
	print "Nb of nodes in the file : %i " % (fixm.getNumberOfNodes())
	fixm.mergeNodes(1e-10)
	print "Nb of non duplicated nodes : %i" % (fixm.getNumberOfNodes())
	# Read and clean Mobile.med
	mobile = ml.MEDFileMesh.New("Mobile.med")
	mobm = mobile.getMeshAtLevel(0)
	mobm.mergeNodes(1e-10)
	# Visualize fixm and mobm with PARAVIEW
	fixm2 = fixm.deepCopy()        # tessellate2D() modifies the current mesh
	fixm2.tessellate2D(0.1)
	fixm2.writeVTK("fixm2.vtu")
	mobm2 = mobm.deepCopy()
	mobm2.tessellate2D(0.1)
	mobm2.writeVTK("mobm2.vtu")
	# mobm2 is in several pieces, take the first one
	zonesInMobm = mobm.partitionBySpreadZone()
	print "Nb of zones in mobm : %i" % (len(zonesInMobm))
	zone1Mobm = mobm[zonesInMobm[0]]
	zone1Mobm.zipCoords()
	displayVTK(zone1Mobm, "zone1Mobm.vtu")
	# Get cell ids from the fix part in the boudning box of zone1Mobm
	ids2 = fixm.getCellsInBoundingBox(zone1Mobm.getBoundingBox(),1e-10)
	partFixm = fixm[ids2]
	partFixm.zipCoords()
	displayVTK(partFixm,"partFixm.vtu")
	# Intersect partFixm with zone1Mobm
	partFixMob, iPart, iMob = ml.MEDCouplingUMesh.Intersect2DMeshes(partFixm,zone1Mobm,1e-10)
	partFixMob.mergeNodes(1e-10)
	# Get the part of partFixm not included in zone1Mobm using partFixMob
	ids3 = iMob.findIdsEqual(-1)
	partFixmWithoutZone1Mobm = partFixMob[ids3]
	displayVTK(partFixmWithoutZone1Mobm,"partFixmWithoutZone1Mobm.vtu")
	# Check that intersection worked properly 
	# Check #0
	areaPartFixm = partFixm.getMeasureField(ml.ON_CELLS).getArray()
	areaPartFixm.abs()
	areaPartFixMob = partFixMob.getMeasureField(ml.ON_CELLS).getArray()
	areaPartFixMob.abs()
	val1=areaPartFixm.accumulate()[0]
	val2=areaPartFixMob.accumulate()[0]
	print "Check #0 %lf == %lf with precision 1e-8? %s" % (val1,val2,str(abs(val1-val2)<1e-8))
	# Check #1
	areaZone1Mobm = zone1Mobm.getMeasureField(ml.ON_CELLS).getArray()
	areaZone1Mobm.abs()
	val3 = areaZone1Mobm.accumulate()[0]
	ids4 = iMob.findIdsNotEqual(-1)
	areaPartFixMob2 = areaPartFixMob[ids4]
	val4 = areaPartFixMob2.accumulate()[0]
	print "Check #1 %lf == %lf with precision 1e-8 ? %s" % (val3,val4,str(abs(val3-val4)<1e-8))
	# Check #2
	isCheck2OK = True
	for icell in xrange(partFixm.getNumberOfCells()):
	    ids5 = iPart.findIdsEqual(icell)
	    areaOfCells = areaPartFixMob[ids5]
	    areaOfCells.abs()
	    if abs(areaOfCells.accumulate()[0] - areaPartFixm[icell]) > 1e-9:
	        isCheck2OK = False
	        pass
	    pass
	print "Check #2? %s" % (str(isCheck2OK))
	# Indicator field creation
	f = ml.MEDCouplingFieldDouble(ml.ON_CELLS,ml.ONE_TIME)
	m = partFixMob.deepCopy()
	m.tessellate2D(0.1)
	f.setMesh(m)
	arr = ml.DataArrayDouble(partFixMob.getNumberOfCells(),1)
	arr[iMob.findIdsEqual(-1)] = 0.
	arr[iMob.findIdsNotEqual(-1)] = 1.
	f.setArray(arr)
	f.checkConsistencyLight()
	f.setName("Zone")
	ml.MEDCouplingFieldDouble.WriteVTK("Zone.vtu",[f])
	# Other zones
	zonesMobm = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords([mobm[zonesInMobm[0]], mobm[zonesInMobm[1]], mobm[zonesInMobm[5]]])
	zonesMobm.zipCoords()
	partFixMob2,iPart2,iMob2 = ml.MEDCouplingUMesh.Intersect2DMeshes(partFixm,zonesMobm,1e-10)
	partFixMob2.mergeNodes(1e-10)
	f2 = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)
	m2 = partFixMob2.deepCopy()
	m2.tessellate2D(0.1)
	f2.setMesh(m2)
	arr = ml.DataArrayDouble(partFixMob2.getNumberOfCells(),1)
	arr[iMob2.findIdsEqual(-1)]=0.
	st = 0
	end = st + len(zonesInMobm[0])
	arr[iMob2.findIdsInRange(st,end)] = 1.
	st += len(zonesInMobm[0]) ; 
	end = st + len(zonesInMobm[1])
	arr[iMob2.findIdsInRange(st,end)] = 2.
	st += len(zonesInMobm[1])
	end = st + len(zonesInMobm[2])
	arr[iMob2.findIdsInRange(st,end)] = 3.
	f2.setArray(arr)
	f2.checkConsistencyLight()
	f2.setName("Zone2")
	ml.MEDCouplingFieldDouble.WriteVTK("Zone2.vtu",[f2])
