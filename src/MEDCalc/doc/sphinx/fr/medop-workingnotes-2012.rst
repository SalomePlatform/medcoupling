.. meta::
   :keywords: maillage, champ, manipulation
   :author: Guillaume Boulant

.. include:: medcalc-definitions.rst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANNEXE: Note de travail concernant le chantier XMED 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. contents:: Sommaire
   :local:
   :backlinks: none


Analyse preliminaire pour le chantier 2012
==========================================

La figure imposée pour le chantier 2012 est l'intégration du nouveau
module de manipulation de champs dans SALOME 6.6 (objectif CEA), en
préparation de la mise en exploitation dans SALOME 7 (objectif EDF).

L'état actuel est:

* Un module SALOME de nom MED intégrant les bibliothèques MEDCoupling,
  MEDLoader, REMAPPER, mais aussi plusieurs packages logiciels
  aujourd'hui obsolètes ou amener à disparaître pour l'échéance
  SALOME7
* Un module SALOME de nom XMED qui fournit les fonctions graphiques
  pour la manipulation de champs.
* Ce module XMED utilise le module VISU pour les vue de contrôle.

La cible est:

* Un module unique (nom à définir, par exemple MEDOP) débarrassé des
  packages logiciels obsolètes et intégrant les fonctions graphiques
  (GUI et TUI).
* L'utilisation du module PARAVIS (au lieu de VISU) pour les vues de
  contrôle.
* L'intégration de MEDCoupling avec YACS (port MED dans YACS par
  exemple).

A examiner:

* voir les attendus concernant les ports MED dans YACS
* interface PARAVIS: utilisation du viewer (et de l'API python) sans chargement du GUI

Tâches de développement
=======================

20120904: Migrer XMED dans MED
------------------------------

Plan de travail:

* Migration des composants + test



20120904: Nettoyage de XSALOME
------------------------------

:status: en cours

* Supprimer les vieilleries de XSALOME:

  - StdHelper -> Basic_Utils (KERNEL)

20120829: mise en place du chantier 2012
----------------------------------------

:status: terminé

L'objectif de cette première étape est de reverser le prototype 2011
(module XMED indépendant) dans la branche V6_main du module MED. On
peut procéder de la manière suivante:

* update de XMED (et XSALOME utilisé par XMED) pour fonctionnement sur
  V6_main
* Eliminer la dépendance à XSALOME
* Supprimer la gestion des multiversion SALOME5/6 au niveau de l'engine

.. warning:: TODO: refaire le point sur les tâches initiées en 2011

