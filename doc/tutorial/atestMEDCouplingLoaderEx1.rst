
.. _python_testmedcouplingloaderex1_solution:

Agitateur - Swirler
~~~~~~~~~~~~~~~~~~~

::

	import MEDLoader as ml
	import numpy as np
	
	# Get available time steps
	data = ml.MEDFileData("agitateur.med")
	ts = data.getFields()[0].getTimeSteps()
	print ts
	# Get position of the swirler
	fMts = data.getFields()["DISTANCE_INTERFACE_ELEM_BODY_ELEM_DOM"]
	f1ts = fMts[(2,-1)]
	fMc = f1ts.getFieldAtLevel(ml.ON_CELLS,0)
	arr = fMc.getArray()
	arr.getMinMaxPerComponent()      # just to see the field variation range per component
	ids = arr.findIdsInRange(0.,1.)
	f2Mc = fMc[ids]
	# Extract pression field on the swirler
	pressMts = data.getFields()["PRESSION_ELEM_DOM"]
	press1ts = pressMts[(2,-1)]
	pressMc = press1ts.getFieldAtLevel(ml.ON_CELLS,0)
	pressOnAgitateurMc = pressMc[ids]
	#
	pressOnAgitateurMc.getMesh().zipCoords()
	# Compute pressure on skin
	agitateurMesh3DMc = pressOnAgitateurMc.getMesh()
	m3DSurf,desc,descI,revDesc,revDescI = agitateurMesh3DMc.buildDescendingConnectivity()
	nbOf3DCellSharing = revDescI.deltaShiftIndex()
	ids2 = nbOf3DCellSharing.findIdsEqual(1)
	agitateurSkinMc = m3DSurf[ids2]
	offsetsOfTupleIdsInField = revDescI[ids2]
	tupleIdsInField = revDesc[offsetsOfTupleIdsInField]
	pressOnSkinAgitateurMc = pressOnAgitateurMc[tupleIdsInField]
	pressOnSkinAgitateurMc.setMesh(agitateurSkinMc)
	# Force field computation
	pressSkin = pressOnSkinAgitateurMc.getArray()
	pressSkin *= 1e5                   # conversion from bar to Pa
	areaSkin = agitateurSkinMc.getMeasureField(True).getArray()
	forceSkin = pressSkin*areaSkin
	normalSkin = agitateurSkinMc.buildOrthogonalField().getArray()
	forceVectSkin = forceSkin*normalSkin
	# Torque computation
	singlePolyhedron = agitateurMesh3DMc.buildSpreadZonesWithPoly()
	singlePolyhedron.orientCorrectlyPolyhedrons()
	centerOfMass = singlePolyhedron.computeCellCenterOfMass()

	barySkin=agitateurSkinMc.computeCellCenterOfMass()
	posSkin = barySkin-centerOfMass

	torquePerCellOnSkin = ml.DataArrayDouble.CrossProduct(posSkin,forceVectSkin)

	zeTorque = torquePerCellOnSkin.accumulate()
	print "couple = %r N.m" % zeTorque[2]
	# Power computation
	speedMts = data.getFields()["VITESSE_ELEM_DOM"]
	speed1ts = speedMts[(2,-1)]
	speedMc = speed1ts.getFieldAtLevel(ml.ON_CELLS,0)
	speedOnSkin = speedMc.getArray()[tupleIdsInField]
	powerSkin = ml.DataArrayDouble.Dot(forceVectSkin,speedOnSkin)
	power = powerSkin.accumulate()[0]
	print "power = %r W"%(power)
	# Eigen vector computation
	x2 = posSkin[:,0]*posSkin[:,0]
	x2 = x2.accumulate()[0]
	y2 = posSkin[:,1]*posSkin[:,1]
	y2 = y2.accumulate()[0]
	xy = posSkin[:,0]*posSkin[:,1]
	xy = xy.accumulate()[0]
	inertiaSkin = np.matrix([[x2,xy],[xy,y2]])
	inertiaSkinValues, inertiaSkinVects = np.linalg.eig(inertiaSkin)
	pos = max(enumerate(inertiaSkinValues), key=lambda x: x[1])[0]
	vect0 = inertiaSkinVects[pos].tolist()[0]
	print vect0

	def computeAngle(locAgitateur1ts):
		fMc = locAgitateur1ts.getFieldAtLevel(ml.ON_CELLS,0)
		arr = fMc.getArray()
		ids = arr.findIdsInRange(0.,1.)
		f2Mc = fMc[ids]
		m3DSurf,desc,descI,revDesc,revDescI = f2Mc.getMesh().buildDescendingConnectivity()
		nbOf3DCellSharing = revDescI.deltaShiftIndex()
		ids2 = nbOf3DCellSharing.findIdsEqual(1)
		agitateurSkinMc = m3DSurf[ids2]
		#
		singlePolyhedron = agitateurMesh3DMc.buildSpreadZonesWithPoly()
		singlePolyhedron.orientCorrectlyPolyhedrons()
		centerOfMass = singlePolyhedron.computeCellCenterOfMass()
		bary = agitateurSkinMc.computeCellCenterOfMass()
		posSkin = bary-centerOfMass
		x2=posSkin[:,0]*posSkin[:,0] ; x2=x2.accumulate()[0]
		y2=posSkin[:,1]*posSkin[:,1] ; y2=y2.accumulate()[0]
		xy=posSkin[:,0]*posSkin[:,1] ; xy=xy.accumulate()[0]
		inertiaSkin = np.matrix([[x2,xy],[xy,y2]])
		inertiaSkinValues,inertiaSkinVects = np.linalg.eig(inertiaSkin)
		pos = max(enumerate(inertiaSkinValues), key=lambda x: x[1])[0]
		vect0 = inertiaSkinVects[pos].tolist()[0]
		return vect0

	vects = len(ts)*[None]
	for itts,locAgitateur1ts in zip(ts,data.getFields()["DISTANCE_INTERFACE_ELEM_BODY_ELEM_DOM"]):
		angle = computeAngle(locAgitateur1ts)
		vects[itts[0]] = angle
		pass

	from math import acos, sqrt
	angle2 = len(ts)*[0.]
	for pos in xrange(2,len(vects)):
	    norm1 = sqrt(vects[pos-1][0]*vects[pos-1][0]+vects[pos-1][1]*vects[pos-1][1])
	    norm2 = sqrt(vects[pos][0]*vects[pos][0]+vects[pos][1]*vects[pos][1])
	    crs = vects[pos-1][0]*vects[pos][0]+vects[pos-1][1]*vects[pos][1]
	    crs /= norm1 ; crs /= norm2 ; crs = min(crs,1.)
	    angle2[pos] = acos(crs) #/(ts[pos][2]-ts[pos-1][2])
	    pass

	omega=sum(angle2)/(ts[-1][2]-ts[0][2])
	print sum(angle2)
	
	print "At timestep (%d,%d) (physical time=%r s) the torque is: %r N.m, power/omega=%r N.m " % (ts[2][0],ts[2][1],ts[2][2],zeTorque[2],power/omega)
