.. meta::
   :keywords: maillage, champ, manipulation, XMED
   :author: Guillaume Boulant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Démonstrateur XMED, vue d'ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Le module XMED est un espace d'expérimentation pour le développement
des opérations de manipulation de champs. Il complète des
développements intégrés directement dans le module MED et gérés dans
la branche CVS BR_medop.

Une maquette est au point pour illustrer les propositions en matière
d'ergonomie d'utilisation et en matière d'architecture technique. La
maquette permet de réaliser des cas d'utilisation de la forme:

* Chargement d'un fichier med dans le module MED (ou publication par
  un code de calcul).
* Sélection graphique des champs de l'étude à mettre à disposition
  dans la console utilisateur ("calculette" en mode texte qui
  concraitement correspond à l'interface python de SALOME).
* Dans la calculette, exécution d'opérations algébriques (+,-,*,/)
  entre champs avec possibilité d'utiliser des scalaires dans des
  opérations de type transformation linéaire (y=ax+b ou y et x sont
  des champs et a et b des scalaires). Egalement quelques fonctions
  mathématiques standard applicables sur des champs (pow, sqrt).
* Possibilité de visualiser les champs produits avec VISU
* Possibilité d'exporter des champs produits dans un fichier med

La figure ci-dessous montre le résultat d'une séquence d'utilisation
dans laquelle les champs "testfield1" et "testfield2" ont été
sélectionnés dans l'arbre d'étude pour être utilisés dans la console
textuelle sous les noms de variables f1 et f2. L'image montre le
contrôle visuel du résultat de l'opération f1+f2-(f1-f2)^2 tapée en
ligne de commande:

.. image:: images/medop-gui-result.png
   :align: center

La séquence ci-après montre le cas d'utilisation complet en
images:

1. Sélection d'un champs sur un pas de temps dans l'arbre d'étude
2. Saisie d'un nom de variable (alias) pour manipuler ce champ. Par
   défaut, le nom du champ est proposé (``testfield1`` ici). Dans
   l'exemple, l'utilisateur remplace par l'alias ``f1``.
3. Contrôle visuel du champ ``testfield1`` manipulé par sa variable
   ``f1`` au moyen de la commande ``f1.visu()``
4. Chargement du champ ``testfield2`` sous le nom ``f2``, exécution de
   l'opération ``f1+f2-(f1-f2)^2`` et contrôle visuel du résultat,
   récupéré ici dans une variable de nom ``result``.

.. |IMG_SELECT| image:: images/medop-gui-selectfield_scale.png
.. |IMG_ALIAS| image:: images/medop-gui-aliasfield_scale.png
.. |IMG_VISU| image:: images/medop-gui-visufield_scale.png
.. |IMG_RESULT| image:: images/medop-gui-result_scale.png

+---------------+---------------+
| |IMG_SELECT|  | |IMG_ALIAS|   |
+---------------+---------------+
| |IMG_VISU|    | |IMG_RESULT|  |
+---------------+---------------+

La solution technique est construite sur les principes suivants:

* Les données MEDMEM sont physiquement chargées par le composant MED,
  c'est-à-dire dans le processus ``Container`` de SALOME, et sont
  référencées dans l'étude SALOME.
* Les opérations sont physiquement des opérations entre objets MEDMEM
  purs qui ont lieu dans le composant MED.
* Les opérations sont pilotées par des objets proxy python instanciés
  dans la console TUI puis manipulés par l'utilisateur. Ces objets
  proxy savent accéder aux objets MEDMEM au travers de leur interface
  CORBA.

Ainsi, l'architecture technique est construite pour pouvoir travailler
sur des données MEDMEM pur en partant de pointeurs CORBA manoeuvrés
depuis des objets python dans l'interface textuelle de
SALOME. L'effort principal a donc porté sur la mise au point de
l'interface technique qui permet de lier des variables représentant
les champs au niveau du GUI (techniquement, la calculette est
l'interpréteur python embarqué dans le GUI, étendu de quelques
fonctions pour la manipulation de champs), alors que les données
MEDMEM sont physiquement disponibles uniquement au niveau des
composants CORBA (et les opérations implémentées dans MEDMEM
uniquement).

Pour le moment, la maquette est limitée à des operations entre champs
qui partagent le même support med (contrainte de MEDMEM) et le
résultat est calculé sur toutes les composantes et tout le domaine de
définition du champs (cette deuxième contrainte est juste parce que
les extentions n'ont pas encore été examinées). Enfin, le support de
gestion des données est supposé être l'étude SALOME et la structure
MED qui y est publiée.
