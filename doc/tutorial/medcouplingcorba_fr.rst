
Visualiser une instance de MEDCoupling dans ParaViS à travers CORBA
-------------------------------------------------------------------

Il peut être intéressant de visualiser directement un maillage ou un champ en mémoire dans 
un process Python avec une session de ParaViS. Cela évite d'avoir à écrire le maillage (ou le champ) sur disque.
Cette technique peut également être utilisée pour :

* faire des noeuds de visualisation dans le module YACS
* maquetter un script Python, et profiter de l'interpreteur interactif Python tout en 
  bénéficiant de l'interface graphique de ParaViS.
  
Nous allons pour ce faire bénéficier des mécanismes de distribution mis en oeuvre dans SALOME sur 
la base du standard `CORBA <http://fr.wikipedia.org/wiki/Common_Object_Request_Broker_Architecture>`_. 
SALOME utilise l'implémentation `omniORB <http://omniorb.sourceforge.net/>`_ et 
`omniORBPy <http://omniorb.sourceforge.net/>`_ du standard.

Début de l'implémentation
~~~~~~~~~~~~~~~~~~~~~~~~~

Pour commencer l'exercice importer le module ``MEDCoupling``
et la classe ``MEDCouplingUMeshServant`` du module Python ``MEDCouplingCorba``. ::

	import MEDCoupling as mc
	from MEDCouplingCorba import MEDCouplingUMeshServant

Créer un maillage
~~~~~~~~~~~~~~~~~

Le but ici est de créer un maillage ``m`` non structuré trivial qui sera envoyé par CORBA à ParaViS. ::

	arr = mc.DataArrayDouble(11)
	arr.iota(0)
	m = mc.MEDCouplingCMesh()
	m.setCoords(arr,arr)
	m = m.buildUnstructured()	

.. note:: Le maillage ``m`` est non struturé, mais s'il est cartésien ça marche aussi !

Créer un servant CORBA
~~~~~~~~~~~~~~~~~~~~~~

Nous allons maintenant créer un *servant* CORBA à partir de ``m``, et faire du process Python courant
un *serveur* CORBA. L'objet ``m`` devient ainsi disponible sur le bus CORBA et pourra être interrogé pour
la visualisation par un service distant, typiquement dans notre cas, le module ParaVis.

Invoquer ``MEDCouplingUMeshServant._this()`` sur ``m`` pour en faire une reférence CORBA (``ref_m``).::

	ref_m = MEDCouplingUMeshServant._this(m)

.. note:: Cette ligne ne se contente pas de faire un servant CORBA mais fait du processus courant Python un serveur CORBA.

Récupérer les identifiants pour ParaViS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ce qui suit est général à tout code omniORBpy. Afficher l'IOR ``ior`` de ``ref_m`` (c'est à dire l'identifiant
unique de l'objet sur le bus CORBA) pour pouvoir passer cette chaîne de caractères au plugin ParaMEDCorbaPlugin 
de ParaViS, et ainsi créer une nouvelle source dans ParaViS. ::

	import CORBA
	orb = CORBA.ORB_init()
	ior = orb.object_to_string(ref_m)
	print ior

Puis, via un copier/coller dans l'IHM ParaViS (Menu "Source -> Para MED Corba Plugin Source"), passer l'IOR. 
On voit s'afficher notre maillage.

Utiliser ParaViS en interactif
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Le but ici est juste de voir le principe. Il s'agit d'un point d'entrée pour réaliser des scripts ParaViS plus évolués.
*Tout en laissant actif ParaViS*, importer le module ``pvsimple`` qui fournit l'interface haut niveau de visualisation : ::

	import salome
	salome.salome_init()
	import pvsimple as pvs

.. note:: Le module ``pvsimple`` est, à peu de choses prêt, identique au module ``paraview.simple``. 
	Il est juste adapté à une utilisation au sein de SALOME. Voir la formation PARAVIS à ce sujet.

Une fois cet import réalisé, le script est automatiquement connecté au même serveur de visualisation que 
l'interface graphique de SALOME (si le module PARAVIS est bien actif !). Nous pouvons donc envoyer des commandes
au serveur de visualisation pour demander l'affichage de notre objet CORBA.  :: 

	# We now talk to the PVServer directly
	import pvsimple as pvs
	pvs.Connect(url)
	src1 = pvs.ParaMEDCorbaPluginSource()
	src1.IORCorba = ior      # This is where we need the CORBA reference of the object created
	dr = pvs.Show(src1)

.. note:: Cela correspond exactement à la manipulation précédente faite via l'interface graphique (ajout d'une nouvelle
	source, saisie de l'IOR, etc ...).

Solution
~~~~~~~~

Le script complet doit être exécuté avec le module PARAVIS (ou une fenêtre ParaView) actif dans l'IHM SALOME!

:ref:`python_testMEDCouplingcorba1_solution`
