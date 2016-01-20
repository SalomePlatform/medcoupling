
.. _python_testMEDCoupling2D_solution:

2D polygons meshing
~~~~~~~~~~~~~~~~~~~

::

	from MEDCoupling import *
	from MEDLoader import *

	from math import *

	numberOfNodes = 25
	numberOfCells = 12

	print "1 ********************"
	spaceDimension = 2

	# Coordinates of central polygon
	X = [1.,0.5,-0.5,-1.,-0.5,0.5]
	Y = [0.,sqrt(3.)/2.,sqrt(3.)/2.,0.,-sqrt(3.)/2.,-sqrt(3.)/2.]

	coordinates = []
	# origin
	coordinates.append(0.)
	coordinates.append(0.)

	# central polygon
	for i in range(6):
		coordinates.append(X[i])
		coordinates.append(Y[i])

	# Coordinates of second couron
	xc = 1.5
	yc = - sqrt(3.)/2.
	d = sqrt(xc*xc+yc*yc)
	a30 = pi/6.0
	a60 = pi/3.0

	for i in range(6):
		angle = a30+a60*i
		xtmp = d*cos(angle)
		ytmp = d*sin(angle)
		start = (i-1)%6
		coordinates.append(xtmp+X[(i-1)%6])
		coordinates.append(ytmp+Y[(i-1)%6])
		coordinates.append(xtmp+X[i%6])
		coordinates.append(ytmp+Y[i%6])
		coordinates.append(xtmp+X[(i+1)%6])
		coordinates.append(ytmp+Y[(i+1)%6])

	print "2 ********************"
	# Creation of mesh
	mesh=MEDCouplingUMesh.New()
	mesh.setMeshDimension(2)
	mesh.allocateCells(numberOfCells)
	mesh.setName("MaFleur")

	myCoords=DataArrayDouble.New()
	myCoords.setValues(coordinates,numberOfNodes,2)
	mesh.setCoords(myCoords)

	print "3 ********************"
	# Connectivity of triangular meshing
	connectivity = []
	for i in range(6):
		connectivity.append(0)
		connectivity.append(i%6+1)
		connectivity.append((i+1)%6+1)
	for i in range(6):
		mesh.insertNextCell(NORM_TRI3,3,connectivity[3*i:3*(i+1)])
		pass

	print "4 ********************"
	# Connectivity of hexagons
	connectivity = []
	for i in range(6):
		start = i%6+1
		connectivity.append(start)
		connectivity.append(start+2*(i+3))
		connectivity.append(start+2*(i+3)+1)
		connectivity.append(start+2*(i+3)+2)
		if i==5:
			connectivity.append(7)
		else:
			connectivity.append(start+2*(i+3)+3)
		connectivity.append((i+1)%6+1)
	for i in range(6):
		mesh.insertNextCell(NORM_POLYGON,6,connectivity[6*i:6*(i+1)])
		pass

	print "5 ********************"
	mesh.checkConsistencyLight()

	medFileName = "MEDCoupling_Fleur.med"
	MEDLoader.WriteUMesh(medFileName,mesh,True)
