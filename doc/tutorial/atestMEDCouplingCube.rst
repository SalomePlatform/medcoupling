
.. _python_testMEDCouplingcube_solution:

3D cube meshing
~~~~~~~~~~~~~~~

::

	from MEDCoupling import *
	from MEDLoader import *
	import MEDLoaderDataForTest

	from math import *

	# Definition of environnement variables
	spaceDimension = 3
	N = 4
	nbOfNodes = N*N*N
	nbOfCells = (N-1)*(N-1)*(N-1)
	nbOfCells2D = (N-1)*(N-1)

	print "1 ********************"
	# Initialisation of coordinates
	coordinates = []
	for k in range(N):
		for j in range(N):
			for i in range(N):
				coordinates.append(float(i))
				coordinates.append(float(j))
				coordinates.append(float(k))
				
	print "2 ********************"
	# Creation of meshing : need following initialisations
	# => Definition of the mesh dimension
	# => Definition of number of cells
	# => Definition of name of meshing
	mesh=MEDCouplingUMesh.New()
	mesh.setMeshDimension(3)
	mesh.allocateCells(nbOfCells+nbOfCells2D)
	mesh.setName("3Dcube")

	print "3 ********************"
	# One remark : only one dimension cells by meshing
	# Construction of volumic meshing
	# => Definition of connectivity
	# => Definition of type of cells
	connectivity = []
	for k in range(N-1):
		for j in range(N-1):
			for i in range(N-1):
				inode = N*N*(k+1)+ N*(j+1)+i
				connectivity.append(inode)
				connectivity.append(inode-N)
				connectivity.append(inode-N+1)
				connectivity.append(inode+1)
				connectivity.append(inode-N*N)
				connectivity.append(inode-N*N-N)
				connectivity.append(inode-N*N-N+1)
				connectivity.append(inode-N*N+1)
	print len(connectivity)
	print 8*(nbOfCells)

	print "4 ********************"
	# Adding cells in meshing
	for i in range(nbOfCells):
		mesh.insertNextCell(NORM_HEXA8,8,connectivity[8*i:8*(i+1)])
		pass

	print "5 ********************"
	# Settings of coordinates and verify if it's OK
	myCoords = DataArrayDouble.New()
	myCoords.setValues(coordinates,nbOfNodes,3)
	mesh.setCoords(myCoords)
	mesh.checkConsistencyLight()

	print "6 ********************"
	# Extraction of surfacic meshing
	pt=[0.,0.,0.]
	vec=[0.,0.,1.]
	nodes = mesh.findNodesOnPlane(pt,vec,1e-12)
	mesh2D = mesh.buildFacePartOfMySelfNode(nodes,True)
	#print mesh2D
	mesh2D.setName("3Dcube")
	mesh2D.checkConsistencyLight()

	print "7 ********************"
	# Creation of field : with following definition
	# => Definition of the mesh support
	# => Definition of field name
	# => Definition of field nature
	field = MEDCouplingFieldDouble.New(ON_CELLS)
	field.setMesh(mesh)
	field.setName("field")
	field.setNature(ExtensiveMaximum)

	# Computing and setting field values
	myCoords=DataArrayDouble.New()
	sampleTab=[]
	bar = mesh.computeCellCenterOfMass()
	print bar.getNbOfElems()
	for i in range(nbOfCells):
		x = bar.getIJ(i+1,1)
		y = bar.getIJ(i+1,2)
		z = bar.getIJ(i+1,3)
		d = sqrt(x*x+y*y+z*z)
		sinus = sin(d)
		#f.setValueIJ(i+1,1,sin(d))
		sampleTab.append(sinus)

	myCoords.setValues(sampleTab,nbOfCells,1)
	field.setArray(myCoords)

	fBF = MEDCouplingFieldDouble.New(ON_CELLS)
	fBF.setMesh(mesh2D)
	fBF.setName("fieldBottomFace")
	fBF.setNature(ExtensiveMaximum)
	Cval = 10.
	myCoords2D=DataArrayDouble.New()
	sampleTab=[]
	for i in range(nbOfCells2D):
		sampleTab.append(Cval)
	myCoords2D.setValues(sampleTab,nbOfCells2D,1)
	fBF.setArray(myCoords2D)

	medFileName = "MEDCoupling_cube3D.med"
	# For note : True / False in Write* functions
	# => True : overwriting existing file
	# => False : add in existing file 
	meshes=[mesh2D,mesh]
	MEDLoader.WriteUMeshes(medFileName,meshes,True);
	MEDLoader.WriteField(medFileName,field,False)
	MEDLoader.WriteField(medFileName,fBF,False)


::

	from MEDCoupling import *
	from MEDLoader import *
	import MEDLoaderDataForTest

	from math import *

	spaceDim3D = 3
	MeshDim2D  = 2
	N = 4
	NbCell2D = (N-1)*(N-1)
	NbCell3D = NbCell2D*(N-1)
	NbNode2D = N*N
	NbNode3D = NbNode2D*N

	# Creation of a extruded meshing
	# input : a 2D meshing and a 1D meshing
	# Creation of 2D meshing
	coordinates = []
	for j in range(N):
		for i in range(N):
			coordinates.append(float(i))
			coordinates.append(float(j))
	Connectivities = [0,4,5,1, 1,5,6,2, 2,6,7,3, 4,8,9,5, 5,9,10,6, 6,10,11,7, 8,12,13,9, 9,13,14,10, 10,14,15,11]
	myCoords = DataArrayDouble.New()
	myCoords.setValues(coordinates,NbNode2D,MeshDim2D)

	m1 = MEDCouplingUMesh.New()
	m1.setMeshDimension(MeshDim2D)
	m1.allocateCells(NbCell2D)
	m1.setCoords(myCoords)
	m1.setName("2D_Support")

	for i in range(NbCell2D):
		m1.insertNextCell(NORM_QUAD4,4,Connectivities[4*i:4*(i+1)])
	m1.changeSpaceDimension(3)

	# Creation of 1D meshing
	coords = [ 0.0, 1.0, 2.0, 3.0 ]
	conn   = [ 0,1, 1,2, 2,3 ]
	m2 = MEDCouplingUMesh.New()
	m2.setMeshDimension(1)
	m2.allocateCells(3)
	m2.insertNextCell(NORM_SEG2,2,conn[0:2])
	m2.insertNextCell(NORM_SEG2,2,conn[2:4])
	m2.insertNextCell(NORM_SEG2,2,conn[4:6])
	myCoords1D=DataArrayDouble.New()
	myCoords1D.setValues(coords,4,1)
	m2.setCoords(myCoords1D)
	m2.changeSpaceDimension(3)

	# Construction of extruded meshing
	center = [0.,0.,0.]
	vector = [0.,1.,0.]
	m2.rotate(center,vector,pi/2.)
	m3 = m1.buildExtrudedMesh(m2,0)
	m3.setName("Extrusion")

	# Construction of group : old fashion mode
	part=[1]
	meshGroup=m3.buildPartOfMySelf(part,True);
	meshGroup.setName("meshGroup");

	medFileName = "MEDCoupling_Extrudedcube3D.med"
	MEDLoader.WriteUMeshesPartition(medFileName,"Extrusion",[m3,meshGroup],True)
	

::

	from MEDCoupling import *
	from MEDLoader import *
	import MEDLoaderDataForTest

	from math import *

	spaceDim3D = 3
	MeshDim2D  = 2
	N = 4
	NbCell2D = (N-1)*(N-1)
	NbCell3D = NbCell2D*(N-1)
	NbNode2D = N*N
	NbNode3D = NbNode2D*N

	# Creation of a grid => Structured mesh
	# Need directions definition
	mesh=MEDCouplingCMesh.New()
	coordsX=DataArrayDouble.New()
	arrX=[ 0., 1., 2., 3. ]
	coordsX.setValues(arrX,4,1)
	coordsY=DataArrayDouble.New()
	arrY=[ 0., 1., 2., 3. ]
	coordsY.setValues(arrY,4,1)
	coordsZ=DataArrayDouble.New()
	arrZ=[ 0., 1., 2., 3. ]
	coordsZ.setValues(arrZ,4,1)
	mesh.setCoords(coordsX,coordsY,coordsZ)
	# Passing structured meshing to unstructured
	# necessary to save meshing
	meshU=mesh.buildUnstructured()
	meshU.setName("Grid")

	# Creation of group : fashion mode
	# if ids cells are known, this step is not to be made
	pt=[1]
	m2 = meshU.buildPartOfMySelf(pt,True);
	ret,tabIdCells = meshU.areCellsIncludedIn(m2,0)
	print ret
	print tabIdCells
	# Definition of the name group
	tabIdCells.setName("meshGroup")

	# Passing MEDCoupling to MEDFile
	fmeshU = MEDFileUMesh.New()
	fmeshU.setName("Grid")
	fmeshU.setDescription("IHopeToConvinceLastMEDMEMUsers")
	myCoords = meshU.getCoords()
	print myCoords
	fmeshU.setCoords(myCoords)
	print "**************************"
	fmeshU.setMeshAtLevel(0,meshU)
	print "**************************"
	fmeshU.setGroupsAtLevel(0,[tabIdCells],False)
	print "**************************"

	medFileName = "MEDCoupling_Gridcube3D.med"
	fmeshU.write(medFileName,2)

