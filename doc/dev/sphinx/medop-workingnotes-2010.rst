.. meta::
   :keywords: maillage, champ, manipulation
   :author: Guillaume Boulant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANNEXE: Note de travail concernant le chantier XMED 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. contents:: Sommaire
   :local:
   :backlinks: none

Principes directeurs du développement
=====================================

En matière de développement:

* On ne cherche pas d'emblée à s'inscrire dans la fabrication d'un
  module SALOME diffusable dans la version d'exploitation 2010 (SALOME
  6). La raison est double: (i) on souhaite au moins pour 2010 ne pas
  devoir tenir compte des contraintes de temps SALOME et (ii) le
  produit envisagé fin 2010 est une maquette qui cherche à éprouver
  l'ergonomie générale d'utilisation et en aucun cas on ne peut
  garantir la réalisation d'un module SALOME compatible avec les
  exigences de mise en exploitation.
* On ne cherche pas d'emblée à capturer tous les cas d'application,
  mais à concevoir un développement qui acceptera les extensions de
  périmètres dans des conditions raisonnables. Aussi, les
  fonctionnalités développées seront celles qui sont nécessaires à la
  réalisation des cas d'application de référence;

En matière d'ergonomie:

* L'interface utilisateur de référence (appelé espace de travail dans
  le volet de spécifications fonctionnelles) est l'interpréteur
  python. Les fonctionnalités doivent être pensées pour un usage
  adapté à une interface textuelle (TUI) de ce type.
* La création d'une interface graphique (GUI) peut être envisagée en
  complément et comme un moyen de manipuler graphiquement les
  fonctionnalités développées pour l'interface textuelle et pour aider
  la préparation des variables dans l'interface python.
* Le modèle d'un processus de manipulation de champs est:

  - Préparation du jeu de variables U, V, ... représentant les champs
    à manipuler. C'est à ce stade que l'on résoud la question de
    sélection des données (dans un champ publié dans l'arbre d'étude,
    par un module de calcul ou par chargement d'un fichier med)
  - Utilisation des variables avec une sémantique la plus proche
    possible du modèle conceptuel et des spécifications
    fonctionnelles;
  - Création des variables qui représentent les résultats des
    fonctions de manipulation;
  - Persistence (fichier med), visualisation (SALOME) ou export (vers
    une structure qui peut être directement utilisable en numpy)

Sur le plan technique:

* On souhaite spécifier clairement le conteneur SALOME des fonctions
  de manipulation de champs. Pour discussion:

  - Il apparaît que les modules SALOME MED et VISU contiennent déjà
    des fonctions qui peuvent faire partie des fonctions de
    manipulations de champs (en particulier pour l'exploration des
    structures MED, leur visualisation et la sélection des données à
    manipuler).
  - Dans la mesure où le module MED n'est pas utilisé à ce jour (en
    tout cas pas sous sa forme de module SALOME) et compte-tenu du
    caractère obsolescent du module VISU (amené à être remplacé sur le
    plan fonctionnel  par le module PARAVIS), on pourrait examiner la
    création d'un module dédié à la manipulation des maillages et des
    champs par l'agrégation technique au sein d'un même module des
    fonctions des modules MED et VISU.

Au moins dans un premier temps, on se donne les limites suivantes:

* Une opération ne peut pas combiner des pas de temps différents. Dans
  l'hypothèse où cette limite venait à être levée, on doit spécifier
  le pas de temps de la donnée résultat;
* Le domaine d'application d'une opération pourra être défini
  exclusivement par la donnée d'un maillage ou un groupe d'éléments du
  maillage;
* On ne traite pas le cas des champs qui prennent leurs valeurs aux
  points de gauss ou aux noeuds par élément. Une particularité de ces
  types de support est que le repérage de la position implique deux
  indices (par exemple l'indice de la maille, puis l'indice du point
  de gauss).

Eléments de conception
======================

Plan général
------------

On peut par exemple imaginer une maquette du genre:

* En C++ dans MEDGUI, charger un fichier med et donner une vue de la
  structure des maillages et des champs dans l'arbre d'étude.
* Sélectionner un élément (par exemple un pas de temps d'un champ) et
  le menu contextuel permet d'exporter ce champ dans la console python
  pour manipulation. Pour cela, s'inspirer de la fonction
  ``XCADGUI::OnLoadScript()`` du XCADGUI pour manoeuvrer un objet
  PythonConsole.
* L'élément est marqué comme ayant été exporté, on peut imaginer une
  récupération ultérieure.
* Exporter un deuxième champ cohérent avec le premier (même pas de
  temps et défini sur le même maillage avec le même support, on
  s'arrange pour).
* Dans la console python, faire les opérations sur les champs
* Publication du champ résultat dans l'arbre d'étude pour sauvegarde
  ultérieure. C'est a priori le gros morceau qui consiste à faire un
  objet CORBA MED à partir d'un objet MED standard, en plus défini
  dans la console python (sous forme d'objet python).

Quand ce premier cas d'utilisation est au point, on peut envisager de
le compléter par les opérations suivantes

* exporter le résultat med dans un fichier
* visualiser les champs produits

Plan de développement:

* Faire une maquette en MEDMEM pur d'abord, car quelque soit le choix
  d'architecture, l'opération physique se déroulera en définitif au
  niveau de MEDMEM pur.
* Prévoir une implémentation des opérations sous forme de fonctions
  informatiques, même les opérations algébriques (+,-,*,/). Pour ces
  dernières et dans certaines conditions (quand on manipule
  directement les strutures MEDMEM et non pas les objets CORBA),
  l'utilisation des formes A+B, A-B, ... peuvent être rendues
  possibles. Dans ce cas, voir la possibilité de combiner plusieurs
  opérations algébriques sur une seule ligne: A+B-C*0.3.
* On peut charger la structure MED sous forme d'objet CORBA publiable
  dans l'étude, de sorte d'avoir accés aux méta-données et pouvoir par
  exemple sélectionner les champs d'intérêt. De cet objet CORBA, on ne
  récupère que les informations nécessaires au chargement d'un champs:
  le nom du champs, le nom de son maillage associé, les identifiants
  du pas de temps, au besoin une structure Field non chargée (par
  exemple pour récupérer plus facilement le maillage).
* Un mécanisme (à développer à partir du PyConsole par exemple)
  pourrait alors permettre le chargement des champs sélectionnés dans
  la console python et sous un nom facile à manoeuvrer. Prendre
  inspiration sur XCADGUI::LoadIntoPythonConsole().
* A priori, les données sont physiquement chargée dans le GUI. Au
  besoin, il semble possible (cf. MED_i::init) de fabriquer une objet
  CORBA field à partir d'un field standard (à tester).

Une autre idée est de récupérer le pointeur CORBA MED dans la console
python et de tirer les données à partir de là. Ajouter une couche de
wrapping python pur pour gérer les cas de simplification (surcharge
des opérations arithmétiques par exemple).

Besoins complémentaires:

* L'interpréteur doit contenir des éléments d'aide (par exemple un
  help qui liste les opérations possibles sur les champs chargés)
* prévoir quelques fonctions de visu et de persistence. Cela commence
  probablement par des fonctions de publication dans l'étude des
  champs créés par les opérations de manipulation. Les champs sont
  physiquement ajouté automatiquement à la structure med par le MedOp
  mais il n'est pas obligatoirement publié => fournir un moyen de
  publication.

Limitations actuelles (liées à la conception de MEDMEM):

* les champs doivent être gérés par la même structure MED car ils
  doivent partager le même support.
* les opérations possibles dans MEDMEM sont entre champs pris sur un
  pas de temps (Q: les pas de temps peuvent-ils être différents).


Développements
--------------

Développement de classes proxy:

* FieldProxy, FieldTimeSeriesProxy
* Attention pour les éries temporelles, le SUPPORT med peut être
  différent en chaque pas de temps (par exemple en cas d'extension
  spatiale du champ au cours du temps).

MEDMEM_MedDataManager:

* FIX: test de l'implémentation C++ au travers de la fonction test() du
  MedOperator ==> OK. Quand on fait la même opération depuis python
  via l'interface SWIG ==> au deuxième appel de getFieldDouble, le
  destructeur du champ semble être appelé. Pb de gestion des pointeurs?


Evolutions à prévoir
====================

Concernant MEDMEM:

* FIX: SALOME_MED::MED::getField devrait pouvoir être appelée
  plusieurs fois de suite puisqu'on recycle la référence si elle est
  déjà chargée.
* IMP: MEDMEM::MED faire une gestion des chargements des champs (par
  exemple avec un getField qui renvoie le champ s'il est déjà chargé
  ou le charge et le renvoie sinon).
* IMP: Récupérer le nom du fichier med à partir de l'objet MED, en
  passant a priori par le driver associé. Plusieurs driver peuvent
  être associés à une structure MED car les données peuvent être
  chargées en plusieurs fois et de plusieurs fichiers. Il faut donc
  étendre la structure MED pour avoir accés à la liste des driver puis
  de cette liste déduire les noms des fichiers.
* IMP: Opérations combinant des champs sur des support différents ne
  peuvent pas être faites par l'API (une exception est levée en cas de
  supports incompatibles), mais on peut imaginer le faire en
  manoeuvrant les tableaux de données directement.
* INF: faire le point sur les fonctions utilitaires autour de MEDMEM
  et de son interface SWIG (ex: dumpMEDMEM.py, med_opfield_test.py).
* IMP: dans MEDMEM::MED et SALOME_MED::MED, pouvoir enlever un champ
  préalablement ajouté: une fonction removeField en complément de
  addField.

Concernant l'interface SALOME_MED:

* IMP: Fonctions algébriques, qui seront implémentées au niveau de la
  structure MED et requêtées au niveau des classes proxy en spécifiant
  les identifiants des champs impliqués et les paramétres requis (pas
  de temps en particulier).

Concernant le module MED:

* IMP: pourvoir exporter la structure med dans un fichier med (la
  structure ayant pu être enrichie par la publication de champs créés
  par les operations de champs.


Historique des travaux
======================

20100726 : mise au point du schéma de conception
------------------------------------------------

Choix entre MEDMEM et MEDCoupling: on reste sur MEDMEM pour plusieurs
raisons:

* MED Coupling ne peut pas gérer des mailles de dimensions différentes
  dans un même modèle (choix faits dans un soucis de performance dans
  l'accès à une structure de donnée compact). On peut contourner le
  problème en définissant deux champs pour traiter chacun des type de
  mailles.
* Un champ repose sur un maillage complet (pas de notion de profil,
  mais cela peut être émulé en créant deux maillages)
* Le concept de point de gauss n'existe pas (pas implémenté)

TODO:

* Idéalement, il conviendrait de faire un état des lieux du module
  MED, en particulier des éléments MEDMEM (le coeur), les interfaces
  CORBA associées (MED.idl implémenté dans le package source
  MEDMEM_I), l'engine (composant SALOME d'interface MED_Gen.idl et
  implémenté dans le package source MED) et le GUI (MedGUI.cxx
  implémenté dans le package source MEDGUI).

* Ergonomie TUI et modèle CORBA associé:

  1. Charger un objet medmem (puis les objets métier mesh et field)
     sur un domaine d'application donné.
  2. En faire des variables disponibles dans l'interface TUI et que
     l'on peut manipuler dans des opérations algébriques.
  3. Pouvoir au besoin en faire des objets CORBA pour l'interface avec
     les autres modules SALOME.

* Compléter le diagramme de la structure informatique de MED (en
  particulier l'implémentation des interface IDL).
* Préparer un module de travail XMED (organisation d'une bibliothèque)

Tests à réaliser:

* Est-il possible de faire des opérations algébriques à partir des
  objets SALOMEMED (objects CORBA MED)?
* Création d'un objet MED_i à partir d'une objet MED pur préalablement
  chargé en mémoire.

A retenir:

* Des opérations de champs sont possibles sur des champs à des pas de
  temps fixés. Si l'opération doit être menée sur plusieurs pas de
  temps, alors itérer sur chaque pas de temps. L'idée ici est
  d'introduire le concept de série temporelle de champs en temps
  qu'objet manipulable.
* Pour deux champs différents de la même structure MED, la données des
  identifiants dt et it ne correspond pas forcément au même instant
  absolu (en tout cas rien ne le garanti, même si c'est tout de même
  une pratique courante).

20101005 : première maquette de démonstration de l'ergonomie en MEDMEM pur
--------------------------------------------------------------------------

XMED: svn révision 16
Travailler avec le fichier de donnée testfield.med joint.


20101007 : Vers une maquette CORBA
----------------------------------

Le contexte d'utilisation des opérations de champs est l'environnement
SALOME. Le support de gestion des données est donc l'étude SALOME. Au
plus bas niveau, les champs sont des objets MEDMEM instanciés dans une
session SALOME (soit par un code de calcul intégré, soit par
chargement des données à partir d'un fichier med). Ces objets sont en
général référencés dans l'étude SALOME sous la forme d'objets CORBA de
classe SALOMEMED::FIELD. Plus exactement, l'étude SALOME gère des
SObject (Study Object) dont un attribut est une référence vers un
objet CORBA de classe SALOMEMED::FIELD qui lui-même encapsule un objet
MEDMEM::Field.

On peut donc envisager une solution dans laquelle on donne à
l'utilisateur des poignées de manipulation des objets
SALOMEMED::FIELD, par exemple au moyen d'un modèle informatique de
type proxy. Cela signifie que l'utilisateur ne manipule pas
directement des objets MEDMEM mais des objets python qui font
l'interface (à concevoir et implémenter, a priori avec un design
pattern de type proxy).

L'utilisation directe des objets MEDMEM aurait pu être une solution
extremement pratique dans la mesure où ces objets en l'état peuvent
être combinés dans des opérations de champs (c'est déjà
implémenté). Par contre, ce procédé souffre de limitations importantes
dans la gestion et la circulation des données pour les différents cas
d'utilisation envisagés (visualisation, export, transfert à un autre
module SALOME).

L'avantage de la solution proposée est multiple:

* Elle permet de travailler sur une structure MED cohérente pour
  intégrer les résultats des opérations de calculs et combiner des
  champs cohérents entre eux. Tout passe par des classes proxy qui
  pourront s'assurer de la cohérence des opérations demandées et
  exécuter automatiquement les fonctions de pré-traitement ou
  post-traitement requises pour ces opérations. On peut imaginer par
  exemple que les requêtes d'opération soient envoyées par les classes
  proxy à la structure MED à laquelle les champs sont associés pour
  piloter l'opération en MEDMEM pur.
* Elle permet d'automatiser un certain nombre d'opérations
  implicites. Par exemple si deux champs ne sont pas définis dans la
  même unité, un changement d'unité peut être effectué automatiquement
  par la classe proxy avant de commander l'opération au niveau
  MEDMEM.
* Elle permet de laisser les données sur le container SALOME et de
  réaliser des opérations sans rappatrier les données en local (qui
  peuvent être en trés grand nombre).
* Elle permet d'étendre facilement l'ergonomie de manipulation des
  champs, par exemple en définissant la notion de *série temporelle de
  champs*, ou encore les concepts de *domaine de définition* évoqués
  dans les spécifications fonctionnelles.
* Elle rend immédiat la circulation des données entre modules SALOME,
  puisque les champs restent accessble par des objets CORBA, en
  particulier pour la visualisation ou l'export des champs produits
  par les opérations.

Elle a cependant des inconvénients et/ou limitations:

* Elle nécessite l'implémentation d'une classe proxy pour encapsuler tous
  les appels aux objets SALOME_MED (et donc MEDMEM). Cette interface
  se limite a priori aux opérations de champs (les opérations
  algébriques dans un premier temps).
* Les champs à manipuler dans une opération donnée doivent être gérés
  par la même structure MED.

Il est à noter également que les interfaces de programmation de
SALOMEMED (interface CORBA pour MEDMEM) devront être étendues pour
permettre des requêtes de manipulations de champs (fonctions addition,
soustraction, multiplication, ...). Pas de contrainte ici sur
l'ergonomie puisque la manipulation par l'utilisateur se fera au
niveau des classes proxy uniquement.


Hypothèses:

* On tente ici une maquette qui exploite dans la mesure du possible le
  fonctionnement actuel du module MED, en particulier la gestion des
  données dans l'étude.
* Dans une deuxième version, on pourra examiner sérieusement la
  révision de la gestion des données dans le module, quitte à la
  spécifier et maquetter dans XMED pour intégration ultérieure dans
  MED. Exemple:

  - Pouvoir gérer plusieurs structures med dans l'étude.

* Enfin, on exploite MEDMEM en l'état. Pour les besoins de la gestion
  des données (gestion des chargements des champs en particulier,
  références croisées pour retrouver le med à partir du champ par
  exemple, ...), il pourra être nécessaire de faire évoluer MEDMEM. Il
  faut pouvoir par ailleurs gérer indifféremment une structure med (et
  les champs qui y sont associés) qu'elle soit créée en mémoire from
  scratch ou chargée d'un fichier (donc attention avec les opérations
  de lecture read(), sur les maillages comme sur les champs). La
  structure med permet d'obtenir les méta données (meta-field par
  exemple) mais ne permet pas de savoir si les données sont
  physiquement chargées ou pas.


Révisions:

* XMED svn revision 21 + tarball MED_SRC-20101014-15h26m.tgz.
  Première version qui permet d'importer un champ dans la console
  python sous la forme d'un FieldProxy. Ne permet pas encore de faire
  des opérations. Introduction dans le module MED de l'interface MEDOP
  pour prendre en charge les opérations sur les champs.


20101019 : Maquette de démonstration pour l'addition
----------------------------------------------------

Cette maquette implémente une solution technique de bout en bout (de
l'interface python aux objets MEDMEM, en passant par le fieldproxy
puis les servants CORBA pour les operations, ...) mais sur le
périmètre de l'addition de champs sur tout leur domaine de définition
et pour un pas de temps donné.

Limitations:

* gére l'addition de champs de type double uniquement (parceque le
  reste n'est pas implémenté)

Révisions:

* XMED: svn révision 25
* MED: cvs tag BR_medop_20101019


20101020: Fonctions complémentaires
-----------------------------------

Cette version test la faisabilité des fonctions complémentaires pour
accompagner la manipulation de champs. Cela comprend en particulier:

* **la sauvegarde des champs produits** dans un fichier med (un champ ou
  toute la structure med). Pour cela, on définit un med proxy comme
  l'extention du SALOME_MED::MED (prévir plutôt d'implémenter ce type
  de fonction au niveau C++ pour permettre un usage au niveau du GUI
  C++?).
* **la visualisation d'un champ** au moyen du module VISU.
* **des fonctions d'aide interactives** pour assister l'utilisateur
  dans la console de manipulation des champs.


Questions:

* peut-on sauvegarder un champ unique?
* peut-on faire en sorte que ce soit l'affectation à une variable qui
  provoque l'ajout du champ à la structure med (ou plus exactement qui
  supprime tous les champs intermédiaires).


Révision:

* XMED: svn revision 31
* MED: cvs tag BR_medop_20101025


20110606: commit avant transfert dans git
-----------------------------------------

* XMED: svn revision 53

Les parties de MED utiles à MEDOP seront reversées dans XMED
dans une première étape, puis le tout dans MED 6 au final. 
