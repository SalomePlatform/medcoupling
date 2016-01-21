
Lecture, écriture d'un fichier MED grâce à l'API avancée de MEDLoader
---------------------------------------------------------------------

L'API avancée de MEDLoader est représentée par les classes ``MEDFile*`` de la bibliothèque MEDLoader.

* Au plus haut niveau, pour l'ensemble du fichier: ``MEDFileData``,
* Pour l'ensemble des maillages du fichier : ``MEDFileMeshes``,
* Pour chacun des maillages : ``MEDFileMeshMultiTS``, ``MEDFileMesh``, ``MEDFileUMesh``, ``MEDFileCMesh``,  
* Pour l'ensemble des champs du fichier : ``MEDFileFields``, ``MEDFileFieldGlobs``, 
* Et enfin pour chacun des champs : ``MEDFileField1TS``, ``MEDFileFieldMultiTS``


Objectif
~~~~~~~~

Ecrire un maillage et un champ à partir de rien, les relire et comparer les résultats.

Points abordés : en utilisant l'API avancée de MEDLoader,

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

Nous créons tout d'abord le même maillage ``targetMesh`` que pour l'API simple. ::

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

Nous construisons ensuite ``targetMesh1`` représentant les sous-constituants (*faces*) du maillage
``targetMesh``, et nous en extrayons seulement les cellules (donc ici des surfaces) [3,4,7,8]. 
Pour plus de détails sur la connectivité descendante, 
consulter la section :ref:`exo-umesh-desc-connec` du deuxième exercise.
Cet ensemble peut par exemple représenter un ensemble d'intérêt pour un calcul : ::

	targetMeshConsti, _, _, _, _ = targetMesh.buildDescendingConnectivity()
	targetMesh1 = targetMeshConsti[[3,4,7,8]]
	targetMesh1.setName(targetMesh.getName())

.. note:: En Python, le underscore ``_`` signifie que l'on attend une valeur de retour, mais qu'on n'en aura pas l'usage 
	(on ne la *bind* pas).
.. note:: ``targetMesh1`` sera sauvé comme étant une partie du même maillage global dans le fichier MED. 
	Il doit donc avoir le même nom. C'est là qu'on voit qu'un maillage au sens MED fichier peut mélanger les dimensions. 

On peut alors écrire les deux maillages dans le fichier "TargetMesh2.med". ::

	meshMEDFile = ml.MEDFileUMesh()
	meshMEDFile.setMeshAtLevel(0,targetMesh)
	meshMEDFile.setMeshAtLevel(-1,targetMesh1)
	meshMEDFile.write("TargetMesh2.med",2)         # 2 stands for 'write from scratch'

Lecture, écriture de groupes de mailles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Créons deux groupes de cellules sur le maillage 2D, c'est à dire au niveau relatif 0 (ici, le niveau relatif 0 correspond
à la 2D, le niveau -1 
correspond à la 1D,  etc ...). Le premier groupe ``grp0_Lev0`` contient les cellules [0,1,3] 
le second ``grp1_Lev0`` les cellules [1,2,3,4] : ::

	grp0_0 = ml.DataArrayInt([0,1,3]) 
	grp0_0.setName("grp0_Lev0")
	grp1_0 = ml.DataArrayInt([1,2,3,4])
	grp1_0.setName("grp1_Lev0")
	meshMEDFile.setGroupsAtLevel(0, [grp0_0,grp1_0])

.. note:: On voit évidemment ici l'importance de nommer les tableaux : c'est le nom qui sera utilisé pour le groupe. 

Créons trois groupes de niveau -1, c'est à dire des groupes de faces. Le premier appelé 
``grp0_LevM1`` aux cellules [0,1], le second appelé ``grp1_LevM1`` aux cellules [0,1,2], et le 3ème ``grp2_LevM1``
aux cellules [1,2,3] : ::

	grp0_M1 = ml.DataArrayInt([0,1])
	grp0_M1.setName("grp0_LevM1")
	grp1_M1 = ml.DataArrayInt([0,1,2])
	grp1_M1.setName("grp1_LevM1")
	grp2_M1 = ml.DataArrayInt([1,2,3])
	grp2_M1.setName("grp2_LevM1")
	meshMEDFile.setGroupsAtLevel(-1,[grp0_M1,grp1_M1,grp2_M1])
	
Ecrivons le tout : ::
	
	meshMEDFile.write("TargetMesh2.med",2)         # 2 stands for 'write from scratch'
	
Nous pouvons ensuite re-lire le fichier MED : ::

	meshMEDFileRead = ml.MEDFileMesh.New("TargetMesh2.med") # a new is needed because it returns a MEDFileUMesh (MEDFileMesh is abstract)
	meshRead0 = meshMEDFileRead.getMeshAtLevel(0)
	meshRead1 = meshMEDFileRead.getMeshAtLevel(-1)
	print "Is level 0 in the file equal to 'targetMesh'?", meshRead0.isEqual(targetMesh,1e-12)
	print "Is level 0 in the file equal to 'targetMesh1'?", meshRead1.isEqual(targetMesh1,1e-12)

Affichons les niveaux disponibles pour le groupe ``grp0_Lev0`` : ::

	print meshMEDFileRead.getGrpNonEmptyLevels("grp0_Lev0")

Et récupérons enfin les identifiants de cellules contenus dans le groupe ``grp0_Lev0`` : ::

	grp0_0_read = meshMEDFileRead.getGroupArr(0,"grp0_Lev0")
	print "Is group 'grp0_Lev0' equal to what is read in the file?" , grp0_0_read.isEqual(grp0_0)

Lire/écrire des champs avec l'API avancée
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Créons un champ de vecteurs simple, aux cellules (P0), avec un seul pas de temps, appelé ``f``. ::

	f = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)
	f.setTime(5.6,7,8)
	f.setArray(targetMesh.computeCellCenterOfMass())
	f.setMesh(targetMesh)
	f.setName("AFieldName")

Stocker ``f`` dans un object ``MEDFileField1TS`` (un champ avec un seul pas de temps -- *one time-step, 1TS*) 
pour préparer l'écriture MED ::

	fMEDFile = ml.MEDFileField1TS()
	fMEDFile.setFieldNoProfileSBT(f)     # No profile desired on the field, Sort By Type

Ajouter le champ au fichier "TargetMesh2.med" ::

	fMEDFile.write("TargetMesh2.med",0) # 0 is paramount to indicate that we *append* (and no overwrite) to the MED file

.. note:: Noter l'utilisation du 0 pour indiquer que nous désirons ajouter au fichier existant.

Lire le champ : ::

	fMEDFileRead = ml.MEDFileField1TS("TargetMesh2.med",f.getName(),7,8)
	fRead1 = fMEDFileRead.getFieldOnMeshAtLevel(ml.ON_CELLS,0,meshMEDFileRead) # Quickest way, not re-reading mesh in the file.
	fRead2 = fMEDFileRead.getFieldAtLevel(ml.ON_CELLS,0)                       # Like above, but this time the mesh is read!
	print "Does the field remain OK with the quick method?", fRead1.isEqual(f,1e-12,1e-12)
	print "Does the field remain OK with the slow method?", fRead2.isEqual(f,1e-12,1e-12)
	
Lire/écrire un champ sur un "profil"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nous allons maintenant voir un concept avancé des fichiers MED, à savoir la possibilité d'écrire un champ sur seulement
une *partie* du maillage. La technique habituellement utilisée est plutôt de mettre des valeurs particulières (e.g. +infinity
soit 1e+300) sur les zones où le champ n'a pas de sens, permettant ainsi de repérer en plus des bugs éventuels lors du calcul.

Le mode de fonctionnement avec les profils reste donc peu courant.

Construisons une réduction aux cellules [1,2,3] de ``f`` et appelons la ``fPart`` : ::

	pfl = ml.DataArrayInt([1,2,3]) 
	pfl.setName("My1stPfl")
	fPart = f.buildSubPart(pfl)
	fPart.setName("fPart")

La stocker dans la structure ``MEDFileField1TS`` et invoquer ``setFieldProfile()``. ::

	fMEDFile2 = ml.MEDFileField1TS()
	fMEDFile2.setFieldProfile(fPart,meshMEDFileRead,0,pfl) # 0 is the relative level (here 0 means 2D)
	fMEDFile2.write("TargetMesh2.med",0) # 0 is paramount to indicate that we *append* (and no overwrite) to the MED file

Lire le champ ``fPart`` du fichier "TargetMesh2.med" et les identifiants de cellules correspondant. ::

	fMEDFileRead2 = ml.MEDFileField1TS("TargetMesh2.med",fPart.getName(),7,8)
	fPartRead, pflRead = fMEDFileRead2.getFieldWithProfile(ml.ON_CELLS,0,meshMEDFileRead)
	print "Is the partial field correclty read?", fPartRead.isEqualWithoutConsideringStr(fPart.getArray(),1e-12)
	print "Is the list of cell identifiers matching?", pflRead.isEqualWithoutConsideringStr(pfl)

Solution
~~~~~~~~

:ref:`python_testMEDLoaderAdvancedAPI1_solution`

