
.. _python_testMEDCouplingRead_solution:

Read med File
~~~~~~~~~~~~~~~~~~~

::

	from MEDCoupling import *
	from MEDLoader import *


	medFileName = "MEDCoupling_cube3D.med"
	MeshName = "3Dcube"
	FieldName = "field"
	Field2DName = "fieldBottomFace"

	# Retrieving meshes
	mesh3D = MEDLoader.ReadUMeshFromFile(medFileName,MeshName,0)
	mesh2D = MEDLoader.ReadUMeshFromFile(medFileName,MeshName,-1)

	# Retrieving fields
	f = MEDLoader.ReadFieldCell(medFileName,mesh3D.getName(),0,FieldName,-1,-1)
	f2 = MEDLoader.ReadFieldCell(medFileName,mesh2D.getName(),-1,Field2DName,-1,-1)

	# Retrieving Coords Mesh
	Coords3D = mesh3D.getCoords()
	Values = Coords3D.getValuesAsTuple()

	# Retrieving field value on 0 tuple
	pos= Values[0]
	res=f.getValueOn(pos)

	# Verify if value is OK
	bar = mesh3D.computeCellCenterOfMass()
	x=bar.getIJ(1,1)
	y=bar.getIJ(1,2)
	z=bar.getIJ(1,3)

	from math import *
	d = sqrt(x*x+y*y+z*z)
	sinus = sin(d)

	if abs(res[0]-sinus)<1.e-5:
		print "OK"
	else:
		print "KO"
