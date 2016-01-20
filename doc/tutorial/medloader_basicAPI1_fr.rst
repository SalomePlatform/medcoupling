
Lecture, écriture d'un fichier MED grâce à l'API basique de MEDLoader
---------------------------------------------------------------------

L'API basique de MEDLoader est contenue dans la classe ``MEDLoader``.
Toutes les méthodes de cette classe sont *statiques* (elles ne dépendent pas d'une instance particulière de la
classe), leurs noms commencent par une majuscule. 
L'ensemble des lectures/écritures sont exécutées à chaque appel de méthode et aucun état interne à la classe n'est
enregistré.

Objectif
~~~~~~~~

Ecrire un maillage et un champ à partir de rien, les relire et comparer les résultats.

Points abordés: en utilisant l'API basique de ``MEDLoader``

* Ecrire un fichier
* Lire un fichier

Début d'implémentation
~~~~~~~~~~~~~~~~~~~~~~

Cet exercice repose comme tous les autres sur le language de script Python. On charge 
le module Python ``MEDLoader``.

Pour information, le module ``MEDCoupling`` complet est inclus dans ``MEDLoader``. Pas besoin de l'importer
si ``MEDLoader`` a été chargé. ::

	import MEDLoader as ml

Lecture, écriture d'un maillage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tout d'abord créons un maillage ``targetMesh`` composé de plusieurs types géométriques. ::

	targetCoords = [-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
	targetConn = [0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
	targetMesh = ml.MEDCouplingUMesh("MyMesh",2)
	targetMesh.allocateCells(5)
	targetMesh.insertNextCell(ml.NORM_TRI3,3,targetConn[4:7])
	targetMesh.insertNextCell(ml.NORM_TRI3,3,targetConn[7:10])
	targetMesh.insertNextCell(ml.NORM_QUAD4,4,targetConn[0:4])
	targetMesh.insertNextCell(ml.NORM_QUAD4,4,targetConn[10:14])
	targetMesh.insertNextCell(ml.NORM_QUAD4,4,targetConn[14:18])
	myCoords = ml.DataArrayDouble(targetCoords,9,2)
	myCoords.setInfoOnComponents(["X [km]","YY [mm]"])
	targetMesh.setCoords(myCoords)

.. note:: Le maillage ``targetMesh`` est ordonné par type géométrique.

Le maillage peut alors directement être écrit ... ::

	ml.WriteUMesh("TargetMesh.med",targetMesh,True)  # True means 'from scratch'

... et relu. ::

	meshRead = ml.ReadUMeshFromFile("TargetMesh.med",targetMesh.getName(),0)
	print "Is the read mesh equal to 'targetMesh' ?", meshRead.isEqual(targetMesh,1e-12)

Lire/Ecrire un champ sur un pas de temps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nous créons maintenant un champ de vecteurs ``f`` aux cellules (P0) avec ``targetMesh`` comme support. 
Ce champ correspond par exemple au temps physique 5.6, repéré par l'itération 7 et la sous-itération 8. 
Nous en profitons pour rappeler
que dans les champs MEDCoupling, le temps physique est donné pour information seulement, le stockage et la plupart des
fonctions de l'API se basent sur les deux derniers entiers. ::

	f = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS, ml.ONE_TIME)
	f.setTime(5.6,7,8)                              # Declare the timestep associated to the field 
	f.setArray(targetMesh.computeCellCenterOfMass())
	f.setMesh(targetMesh)
	f.setName("AFieldName")
	ml.WriteField("MyFirstField.med",f,True)

Question subsidiaire : à quoi correspond le champ ainsi créé ?

.. note:: Le maillage **et** le champ sont écrits d'un seul coup dans le fichier "MyFirstField.med".

Nous relisons ensuite MyFirstField.med : ::

	f2 = ml.ReadFieldCell("MyFirstField.med", f.getMesh().getName(), 0, f.getName(), 7, 8)
	print "Is the read field identical to 'f' ?", f2.isEqual(f,1e-12,1e-12)
	
.. note:: Lors de la lecture du champ, on doit donc connaître: son nom, le nom de sa mesh de support
	et le pas de temps voulu. Des fonctions du type ``MEDFileFields.getFieldsNames()`` ou encore 
	``MEDFileMeshes.getMeshesNames()`` aident à cela.
	
.. note:: Le nom ``ReadFieldCell()`` rappelle que le champ doit être lu aux cellules. Souvenez-vous que suivant la 
	norme MED fichier, un même champ peut avoir une partie de ses données stockées aux cellules, mais aussi 
	simultanément aux noeuds, aux points de Gauss, etc ... même si ce genre de mélange exotique n'est généralement
	pas conseillé.

Lire/Ecrire un champ sur plusieurs pas de temps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ici contrairement au cas précédent, nous écrivons en plusieurs fois dans le *même* fichier MED.
Ecrivons tout d'abord le maillage. ::

	ml.WriteUMesh("MySecondField.med",f.getMesh(),True)
	
Ensuite, nous écrivons seulement les informations relatives au champ (principalement son tableau de valeurs en fait
). ::

	ml.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f)   # mesh is not re-written
	
Nous rajoutons ensuite un second pas de temps sur le *même* maillage. ::

	f2 = f.clone(True)         # 'True' means that we need a deep copy  
	f2.getArray()[:] = 2.0
	f2.setTime(7.8,9,10)
	ml.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f2)

Maintenant le fichier "MySecondField.med" contient le maillage et un champ à deux pas de temps porté par ce maillage.

Nous pouvons relire tout cela avec des méthodes similaires à ce qui a été vu précédemment : ::

	f3 = ml.ReadFieldCell("MySecondField.med",f.getMesh().getName(),0,f.getName(),7,8)
	print "Is the field read in file equals to 'f' ?", f.isEqual(f3,1e-12,1e-12)
	f4 = ml.ReadFieldCell("MySecondField.med",f.getMesh().getName(),0,f.getName(),9,10)
	print "Is the field read in file equals to 'f2' ?", f2.isEqual(f4,1e-12,1e-12)

Solution
~~~~~~~~

:ref:`python_testMEDLoaderBasicAPI1_solution`
