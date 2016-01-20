
Reading, Writing a MED file using MEDLoader's basic API
-------------------------------------------------------

The basic API is incarnated by the MEDLoader class.
All methods in that class are static.
The sessions of read/write are done on each call of a method.

Objective
~~~~~~~~~

Write a mesh and a field from scratch, re-read them and compare the result.

Topics covered:
* Write using MEDLoader's basic API
* Read using MEDLoader's basic API

Implementation start
~~~~~~~~~~~~~~~~~~~~

To implement this exercise we use the Python scripting language and import the MEDLoader Python module.
The whole MEDCoupling module is fully included in MEDLoader. No need to import MEDCoupling when MEDLoader has been loaded. ::

	import MEDLoader as ml

Writing/Reading a mesh
~~~~~~~~~~~~~~~~~~~~~~

First of all, creation of a mesh "targetMesh". ::

	targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        targetMesh=ml.MEDCouplingUMesh.New("MyMesh",2)
        targetMesh.allocateCells(5)
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
	targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        myCoords=ml.DataArrayDouble.New(targetCoords,9,2)
	myCoords.setInfoOnComponents(["X [km]","YY [mm]"])
        targetMesh.setCoords(myCoords)
        
.. note:: targetMesh is ordered by geometric type.

We are then ready to write it. ::

	ml.WriteUMesh("TargetMesh.med",targetMesh,True)

Then trying to read it. ::

	meshRead=ml.ReadUMeshFromFile("TargetMesh.med",targetMesh.getName(),0)
	print "Is the mesh read in file equals targetMesh? %s"%(meshRead.isEqual(targetMesh,1e-12))

Writing/Reading a field on one time step at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Creation of a vector field "f" on cell supported by "targetMesh". ::

	f=ml.MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
	f.setTime(5.6,7,8)
	f.setArray(targetMesh.computeCellCenterOfMass())
	f.setMesh(targetMesh)
	f.setName("AFieldName")
	ml.WriteField("MyFirstField.med",f,True)

.. note:: Mesh AND Field is written at once into MyFirstField.

Reading into MyFirstField.med ::

	f2=ml.ReadFieldCell("MyFirstField.med",f.getMesh().getName(),0,f.getName(),7,8)
	print "Is the field read in file equals f ? %s"%(f2.isEqual(f,1e-12,1e-12))

Writing/Reading a field on one or many times steps in "multi-session mode"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here contrary to the previous steps, we are going to write in a multi-session mode on the same MED file.
First dealing with the mesh. ::

	ml.WriteUMesh("MySecondField.med",f.getMesh(),True)
	
Then writing only array part of field. ::

	ml.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f)
	
Then put a another time step. ::

	f2=f.clone(True)
	f2.getArray()[:]=2.0
	f2.setTime(7.8,9,10)
	ml.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f2)

Now "MySecondField.med" file contains 2 time steps.

Solution
~~~~~~~~

:ref:`python_testMEDLoaderBasicAPI1_solution`
