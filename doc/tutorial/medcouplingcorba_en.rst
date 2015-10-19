
Visualize a MEDCoupling instance in ParaViS through CORBA
---------------------------------------------------------

ParaViS can be used to directly visualize a mesh or a field stored in memory in a Python 
process. For information, this technique will become the preferred choice for the MED
Calculator in a future Salome release. 
The following use cases can also be mentioned:

* YACS, to create visualization nodes
* create a Python mock-up script and use the standard Python interpreter whilst benefiting
  from the ParaViS graphical interface

Implementation start
~~~~~~~~~~~~~~~~~~~~

Import the whole Python module MEDCouplingCorba. ::

	from MEDCouplingCorba import *


Create a 2D MEDCouplingUMesh instance 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a trivial unstructured  mesh "m" which will be sent through CORBA to ParaViS.
::

	arr=DataArrayDouble(11)
	arr.iota(0)
	m=MEDCouplingCMesh()
	m.setCoords(arr,arr)
	m=m.buildUnstructured()	

.. note:: "m" is unstructured but a Cartesian mesh would also work perfectly fine.

Create a CORBA servant from "m", and turn the Python process into a CORBA server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Invoke MEDCouplingUMeshServant._this() on "m" to turn it into a CORBA reference ("ref_m").
::

	ref_m=MEDCouplingUMeshServant._this(m)

.. note:: This command doesn't only create a CORBA servant but also makes the current 
	Python process a full CORBA server.

Read the identifiers that are passed to ParaViS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

What follows holds for any omniORBpy code. Display the IOR "ior" of "ref_m".
This character string is given to the ParaViS  plugin (ParaMEDCorbaPlugin) to create 
a new ParaViS source.
::

	import CORBA
	orb=CORBA.ORB_init()
	ior=orb.object_to_string(ref_m)
	print ior

A simple copy/paste in the ParaViS GUI allows to create the source and to have our
mesh rendered on screen.

Use ParaViS interactively
~~~~~~~~~~~~~~~~~~~~~~~~~

This section simply highlights what can be done in principle. It should be regarded
as a starting point towards the creation of more advanced scripts.
With ParaViS still up, retrieve a remote handle on ParaViS:
::

	import PARAVIS_Gen_idl
	import salome
	salome.salome_init()
	paravis=salome.lcc.FindOrLoadComponent("FactoryServer","PARAVIS")

Then send a script to ParaViS so that it displays "m":
::

	script="""
	src1 = ParaMEDCorbaPluginSource()
	src1.IORCorba = '%s'
	asc=GetAnimationScene()
	rw=GetRenderView()
	dr=Show()\ndr.Visibility = 1
	Render()
	"""
	content=script%(ior)
	paravis.ExecuteScript(content)


Solution
~~~~~~~~

:ref:`python_testMEDCouplingcorba1_solution`
