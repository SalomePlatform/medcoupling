
Reading a med file
-------------------

Objective
~~~~~~~~~

The MEDLoader class also allows to read a med file. 

We will use the case of the 3D cube in order to retrieve all meshes and fields along with value on one node.

Beginning of implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To implement this exercice we use the python language script and import the MEDCoupling and MEDLoader parts of the MED module. We need also mathematical functions, so we import the python math module::

	from MEDCoupling import *
	from MEDLoader import *
	from math import *


Then define some variables::

	medFileName = "MEDCoupling_cube3D.med"
	MeshName = "3Dcube"
	FieldName = "field"
	Field2DName = "fieldBottomFace"

Retrieving 3D mesh and associated field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to use MEDLoader API in order to read med file. Read functions need to give the real dimension of the mesh to max.
This information is given by a number : 0,-1 or -2.

 * 0 means  the high dimension of the mesh
 * -1 means the second high dimension of the mesh
 
and the iteration and order of the field. In our case, since there is no iteration, it's -1 for these 2 arguments::

	mesh3D = MEDLoader.ReadUMeshFromFile(medFileName,MeshName,0)
	f = MEDLoader.ReadFieldCell(medFileName,mesh3D.getName(),0,FieldName,-1,-1)


Retrieving 2D mesh and associated field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do the same thing for the 2D mesh and the associated field::

	mesh2D = MEDLoader.ReadUMeshFromFile(...)
	f2 = MEDLoader.ReadFieldCell(...)

Retrieving mesh coords
~~~~~~~~~~~~~~~~~~~~~~

::

	# Retrieving Coords Mesh
	Coords3D = mesh3D.getCoords()
	Values = Coords3D.getValuesAsTuple()

Retrieving field values
~~~~~~~~~~~~~~~~~~~~~~~~

::

	# Retrieving field value on 0 tuple
	pos= Values[...]
	res=f.getValueOn(pos)

	# Verify if value is OK
	bar = mesh3D.computeCellCenterOfMass()
	x=bar.getIJ(...)
	y=bar.getIJ(...)
	z=bar.getIJ(...)

	from math import *
	d = sqrt(x*x+y*y+z*z)
	sinus = sin(d)

	if abs(res[0]-sinus)<1.e-5:
		print "OK"
	else:
		print "KO"

Solution
~~~~~~~~

:ref:`python_testMEDCouplingRead_solution`
