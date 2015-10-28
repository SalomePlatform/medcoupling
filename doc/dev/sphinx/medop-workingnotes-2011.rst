.. meta::
   :keywords: maillage, champ, manipulation
   :author: Guillaume Boulant

.. include:: medcalc-definitions.rst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANNEXE: Note de travail concernant le chantier XMED 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. contents:: Sommaire
   :local:
   :backlinks: none

Cas d'utilisation métier
========================

On illustre par un exemple (Christophe Vallet, R&D/MMC, 1/7/2011)::

 J'ai souvent des fichiers med de résultats de calcul, et j'aimerais y
 ajouter de nouveaux champs issus de champs existants. J'aimerais
 aussi pouvoir créer de nouveaux meds plus petits par extraction de
 certaines composantes de champs, certains groupes ou certains pas de
 temps.

On peut exprimer le besoin sous la forme des cas d'utilisation
suivants (use cases):

* **UC1: combiner dans un même fichier med des champs issus de
  plusieurs sources de données**. On peut par exemple charger un
  premier fichier, puis ajouter à cette base des champs issus d'autre
  fichiers ou générés par manipulation de champs, ou encore générés
  par un module de calcul qui produirait directement du MEDCoupling.
* **UC2: créer un champ contenant certaines composantes d'un autre
  champ**. On pense ici aux fonctions de restriction, qui permettraient
  de récupérer certaines composantes uniquement.
* **UC3: créer un champ contenant certains pas de temps d'un autre
  champ**. C'est un cas particulier des fonctions de restriction
  évoquées ci-dessus.
* **UC4: créer un champ comme la limitation d'un autre champ à un
  groupe de mailles**. C'est un cas particulier des fonctions de
  restriction évoquées ci-dessus. Notion de domaine spatial. A
  priori la notion de groupe est définie dans MEDLoader.

On peut ajouter également les UC identifiés pour la maquette 2010:

* **UC5: comparer des champs issus de source de données différentes**,
  par exemple des champs chargés de deux fichiers med différents et
  qui s'appuient sur le même maillage (au moins conceptuellement).  Le
  problème technique ici est de pouvoir changer le maillage d'un
  champ, pour ramener tous les champs sur le même maillage (au sens
  informatique). Ceci est une contrainte de MEDCoupling, les
  opérations sur des champs A et B imposent que A et B soient définis
  sur le même maillage, i.e. le même objet informatique.
* **UC6: créer un champ de toute pièce sur un maillage**, ou un groupe
  de mailles. Ce cas d'usage est typiquement prévu pour produire les
  conditions de chargement initial d'une structure. Il s'agit ici
  d'initialiser un champ à partir de zéro sur une surface prédéfinie
  de la géométrie (par exemple spécifiée par un nom de groupe de
  mailles).

Pour UC5: les sources de données sont référencées dans l'object
browser. On importe explicitement les données dans l'espace de
travail. On peut détecter que les maillages sont identiques et on
propose à l'utilisateur de transférer le champ sur le maillage déjà
présent. Sinon, les champs devront être référencés sur des maillages
distincts dans l'arbre de l'espace de travail.

Analyses préliminaires pour le chantier 2011
============================================

On fait le choix pour le chantier 2011 de travailler à partir de la
bibliothèque MEDCoupling (et non plus MEDMEM comme c'était le cas dans
le démonstrateur 2011).

Analyse de MEDCoupling et MEDLoader
-----------------------------------

MEDCoupling est l'implémentation du modèle de données MED (avec
recherche de minimisation des dépendances logicielles) et MEDLoader
fournit une ensemble de fonctions pour le chargement des structures
MEDCoupling depuis un fichier ou inversement leur sauvegarde sous
forme de fichiers.

Dans l'implémentation MEDCoupling, un champ est l'ensemble des valeurs
d'une grandeur physique sur un maillage pour un pas de temps donné. Un
champ est caractérisé par:

* un support spatial, le maillage
* un type de discrétisation spatial, défini par l'emplacement des
  valeurs sur le maillage (sur les noeuds, sur les cellules, aux
  points de gauss, ...) et le mode d'interpolation spatial (P0, P1,
  etc)
* un pas de temps, défini par deux entiers (iteration, order) et un
  réel (timestamps)

Dans cette implémentation, il existe une association 1..n entre un
maillage et un champ (alors que dans MEDMEM, la structure
intermédiaire SUPPORT est implémentée).

MEDCouplingCorba fournit un ensemble de servants CORBA pour manoeuvrer
des structures MEDCoupling au travers du bus CORBA. L'interface à ce
jour est délibérément réduite. Des classes dites "Cliente" sont
fournies pour piloter les servants CORBA depuis un contexte
client. Par exemple ``MEDCouplingFieldDoubleClient`` fournit une
fonction de création d'une structure MEDCoupling à partir d'un
pointeur vers un servant CORBA. La structure est créée localement
(dans le contexte client) avec duplication des données issue de la
structure encapsulée par le servant CORBA (récupération par la
fonction de sérialisation).

Aucune interface CORBA n'est défini pour MEDLoader.

Questions:

* Voir comment sont créés les servants, et surtout comment ils sont
  récupérés (via le lcc?)
* Comment peut-on définir un champ sur un groupe de mailles (et non
  pas sur le maillage complet)? Comment peut-on extraire le champs
  circoncit à une groupe de mailles pour des opérations.

  - R: méthode changeUnderlyingMesh

* Comment manipuler deux champs chargées de fichiers différents mais
  construit sur le même maillage (conceptuellement). On peut forcer la
  réassociation d'un champ sur un autre maillage?
* Manipuler des champs de pas de temps différents? Différentes
  composantes d'un ou plusieurs champs?
* Comment importer un MedCoupling dans PARAVIS? (dans VISU?)?

* mapper sur une image

Improvments:

* MEDLoader::Write should raise an exception if the filepath is not writable
* MEDDataManager: développer une classe chapeau sur MEDCoupling et
  MEDLoader pour  aider au chargement et la gestion de données MED
  (orienté manipulation de champs). Cette classe serait associée des
  structures légères FieldHandler et MeshHandler et des listes
  correspondantes pour la navigation dans les méta-données.
* Sur base du MEDDataManager, prévoir des ports med pour yacs par
  lesquels pourrait transiter des handler.

Nouveaux concepts à prendre en compte
-------------------------------------

Au démarrage du chantier 2011, on observe que les concepts suivants
sont introduits dans le module MED:

* Le conteneur MED n'existe plus, utiliser MEDFILEBROWSER pour charger
  les fichiers med et obtenir les informations générales sur le
  contenu.
* MEDFILEBROWSER: remplace le concept de driver et fournit les
  fonctions précédemment fournies par la classe MED pour obtenir les
  informations de structure.
* Concept d'Extractor pour une lecture sélective des données de champs
  (suivant un critère d'extraction)
* Il n'est plus nécessaire d'appeler les méthodes read explicitement
  sur les objets (MESH et FIELD) pour charger les données. Par
  ailleurs, on peut définir deux fois le même champs (double
  chargement a priori) sans lever d'exception).


Analyse de conception pour le chantier 2011
===========================================

Composants SALOME (interfaces IDL)
----------------------------------

* MEDDataManager: défini une structure FIELD pour identifier un champ
  dans les requêtes. Il s'occupe également de la récupération physique
  des données, quelqu'en soit la source (fichier avec MEDLoader, autre
  module SALOME comme PARAVIS avec une méthode à définir)
* MEDCalculator: s'occupe des requêtes de calcul dont les arguments sont
  les structures FIELD du MEDDataManager. Reprendre l'interface de
  MEDOP.

Use case à réaliser depuis un client python:

* UC01: ajouter un fichier d'entrée et accéder aux informations
  concernant les champs. Ex: récupérer une structure champs par la
  donnée des paramètres primaires (nom identifiant, dt, it, nom du
  maillage).
* UC02: créer des champs et les ajouter au MEDDataManager
* UC03: mener des opérations basique sur les champs en console python

Interface Utilisateur
---------------------

L'interface utilisateur est composée des parties suivantes:

* une partie GUI (appelée par la suite MEDGUI) qui s'occupe de piloter
  le chargement des données dans l'espace de travail, au moyen d'une
  interface graphique;
* une partie TUI (appelée par la suite MEDTUI) qui s'occupe de piloter
  la création de champs, au moyen de commandes exécutées dans la
  console python.

Le principe est que les champs sont préalablement chargés au niveau du
composant SALOME au moyen de l'interface graphique (MEDGUI), puis
manoeuvrés depuis l'application SALOME au moyen de variables proxy
définies dans la console python (MEDTUI). Au chargement, les champs
sont indéxés par le MEDDataManager, puis les index sont rendus
accessibles au niveau du GUI au moyen d'une représentation
arborescente de la structure MED. Les feuilles de l'arbre
correspondent à des champs qui peuvent être sélectionnés et dont
l'index peut être obtenu de la sélection.

L'espace de travail est organisé autour du concept de
"workspace". L'étude SALOME liste les datasource (les fichiers source
des données med, mais peut-être aussi les référence vers des objets
MED déjà existants ou chargé dans PARAVIZ). Une vue complémentaire
permet de voir la structure fine d'une source de données.

Concernant MEDGUI:

* la représentation des données (les champs et les maillages associés)
  doit permettre de récupérer par l'interface graphique les
  identifiants des champs à manipuler (a priori les structures FIELD
  définies par le composant MEDDataManager). Cela conduit à la mise en
  place des composants suivants:

  - MedDataModel hérité de TreeData. Il est peuplé avec les
    méta-données décrivant la structure MED explorée.
  - MedGuiManager qui permet l'implantation du doc widget de
    présentation

TODO:

* specifier le concept de workspace (qui a une entrée dans l'étude?)
  en bijection avec un datamanager
* identifier des interlocuteur/utilisateur pour l'aspect ergonomie d'usage

Concernant MEDTUI:

* Il fournit les classes FieldProxy

Questions:

* Comment traiter le cas du travail sur des composantes ciblées, plus
  généralement, comment introduire le concept de domaine
  d'application?
* Prévoir des fonctions génériques (initialisation d'un champ sur un
  maillage avec une fonction analytique de la position, sauvegarder
  les champs créés dans un fichier med)


Tâches de développement
=======================

T20110622.1: Gestion des données internes
-----------------------------------------

**Status: terminé.**
Suite: fonction de sauvegarde au niveau graphique également

On vise les cas d'utiliation suivants:

* UC1: intégrer dans le datamodel du gui un champ créé dans la console
  python (et donc présent dans le datamanager du composant). Définir
  l'utilité?
* UC2: renommer un champ et plus généralement changer ses méta-données
  (avec assurance de synchronisation entre toutes les données).
* UC3: sauvegarder une sélection de champs. La sélection peut se faire
  dans l'arbre du datamodel gui.

WARN: robustesse de fieldproxy



T20110622.2: UC Initialisation/Création de champs
-------------------------------------------------

**Status: à faire**

Les cas implémentés à ce jour sont la création de champs à partir de
champs existants et chargés d'un fichier med. On souhaite ici réaliser
des cas 'utilisation autour de la création de champs "from scratch",
s'appuyant tout de même sur un maillage chargé.

UC01: Sélection d'un groupe de maille dans SMESH pour initialiser un
champ (par exemple les conditions limites d'un problème de calcul).

UC02: créer un champ avec des restrictions qui définissent le domaine
d'application des opération de champs.

UC03: créer un champ à partir d'une image (codes rgb utilisé comme les
composantes du champs vectoriel ou niveaux de gris pour un champ
scalaire. Attention, pour ça, il faudra a priori fiare une projection
du maillage cartesien de l'image sur le maillage (quelconque) sur
lequel on souhaite définir le champ.

UC04: créer un champ à partir d'un tableau numpy

De manière générale, ce type de création sera assisté par le
MEDGUI. Au niveau MEDTUI, les fonctions pourraient être fastidieuses
pour l'utilisateur.

Par exemple, prévoir un menu contextuel qui propose les opérations
possibles en fonction de la sélection (en plus de la fonction d'import
dans la console python).

TODO:

* développer les fonctions d'initialisation, par exemple au moyen
  d'applyFunc et du mécanisme de callable?

T20110622.3: documentation contextuel
-------------------------------------

**Status: à faire**

* Remettre toutes les commandes dans le même fichier (fusionner cmdtools
  et fieldtools)
* Faire un modèle générique de command (classe de base
* Batir la doc des commandes sur cette base (lister toutes les
  instances de type Command par exemple)

T20110622.4: remontée des exception du composant MEDCalculator
--------------------------------------------------------------

**Status: en cours, compléter la couverture**

Pour des messages contextuel sur les erreurs de calcul (ex: division
par 0)

* Poursuivre le travail fait sur getMedEventListener
* Protéger tous les appels au composants effectués depuis la console
  python (prendre example sur la commande save)

T20110624.1: gestion des données GUI
------------------------------------

**Status: à faire**



Le workspace a une entrée dans l'obrowser. Sur cette entrée on peut:

* supprimer: supprime tout les champs associés
* sauvegarder. Dans ce cas, on rappelle l'ensemble des champs pour
  cocher ceux qu'on veut sauvegarder.

Le gui data model est réservé aux opérations sur les champs et à
piloter leur import dans la console python.

TODO:

* Spécifier les concepts de workspace, database, et datasource, espace
  de gestion, ... et les associations. Simplifier avec l'appuie de use
  cases.
* Mécanisme de mise à jour du TreeView de XSALOME (aujourd'hui, seul
  l'ajout addChild est implémenté
* Clic droit sur objets de l'arbre: dans la notification TreeView ->
  WorkspaceController, faire remonter l'évènement clic droit ainsi que la
  liste des éléments sélectionné pour faire générer le menu contextuel
  au niveau du WorkspaceController qui peut déterminer le contexte métier
  (le TreeView ne le connaît pas).
* Définir des DataObject pour les maillages, les séries temporelles et
  les champs


Spécification des espaces de données:

* MEDDataManager dépend de l'étude (pour permettre la publication
  d'information dans une étude SALOME).
* créer "sourcid = MEDDataManager::addDataSource(filename)", suivie de
  requetes getFields(sourceid), getMeshes(sourceid)
* les espaces de données: dataspace, workspace. Un seul workspace par
  étude, mais autand de datasources que l'on souhaite dans le
  dataspace. Les datasources sont rangés dans l'étude (le dataspace)
  et sont non modifiables après chargement (référence des sources de
  données).


T20110628.1: extention à d'autres objets SALOME
-----------------------------------------------

**Status: suspendu**

On doit reposer la question de l'existance de l'arbre indépendant
(DockWidget), d'une part, et l'extention aux autres objets (GEOM et
SMESH en particulier) du principe de sélection graphique pour
utilisation dans la console python, d'autre part.


T20110628.2: visualisation d'un champ avec PARAVIS
--------------------------------------------------

**Status: terminé (pour une première version)**
Suite: de nombreux défauts subsistent

Questions/remarques:

* Pb au démarrage du module: VisTrails fails to start
* Peux-t-on piloter la vue 3D sans charger le module? (voir
  myparavis.py)
* Comment donner un nom au MEDReader1 dans l'arbre Pipeline?
* Comment utiliser directement les objets MEDCouplingField?


T20110706.1: documentation du module
------------------------------------

**Status: en cours (10%)**

Documenter les commandes TUI puis l'utilisation générale de
l'interafce graphique. Mentionner l'existance de la commande medop.sh
pour travailler exclusivement en mode texte (utile pour les tests
rapides).

Documenter les modalités d'exécution des tests.

T20110708.1: helper python pour MEDCoupling
-------------------------------------------

**Status: en attente (pas urgent)**

Faire un helper python dans le package xmed qui permet de faire du
medcoupling facilement (essentiellement pour simplifier le chargement,
puis la sélection des données). Cela demanderait de faire un
MedDataManager comme une class C++ pure (non CORBA). Cette classe
travaillerait par exemple uniquement avec des id et des liste d'id, et
fournirait des fonctions d'affichage (comme le ``ls`` et le ``la``)
pour obtenir des meta-information.

Le servant MedDataManager pourrait être une surcouche de cette classe
c++ pure.

T20110708.2: analyses et tests
------------------------------

TODO:

* créer un fichier de test avec plusieurs pas de temps
* créer un fichier de test avec des groupes de mailles


T20110728.1: refactoring MEDDataManager
---------------------------------------

Refactoring pour une meilleur association entre FieldHandler et MeshHandler:

* dans la mesure du possible utiliser les id plutôt que les handler en
  arguments des fonctions d'appel des objets
* A chaque champ (FieldHandler), on doit associer un meshid (et de
  manière optionnelle un fieldseriesId, si le champ peut être associé
  à une serie temporelle. A priori faisable uniquement au chargement
  du datasource).
* Pour cela, revoir les fonctions internes newFieldHandler et addField
  ou prévoir de les compléter à chaque fois qu'elles sont appelée avec
  les informations concernant le meshid.
* addField est utilisée par le MEDCalculator
* Attention au raffraichissement des données handler au niveau du
  Workspace. Peut-être le mieux est que les fieldproxy contiennent
  uniquement le fieldid, et qu'ils interroge le datamanager à chaque
  fois qu'ils ont besoin d'une donnée. Voir aussi les notifications
  via le MEDEventListener?  **Le plus simple est de faire la mise à
  jour lors de l'appel à la méthode __repr__ du fieldproxy, i.e. quand
  on essaye d'afficher les données**. Parceque sinon il n'y a pas de
  problème puisque que le calculateur travaille à partir des id.


Petites améliorations du DataspaceController:

* Au OnUseInWorkspace, stocker (dans la mesure du possible) le nom de
  l'alias python dans un attribut du sobject.
* Dans DlgChangeUnderLyingMesh, expliquer que le champs sera dupliquer
  est posé dans le WS. On peut donc proposer en option de lui associer
  un alias pour manipulation dans la console



