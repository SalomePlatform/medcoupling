.. meta::
   :keywords: SALOME, development
   :author: Guillaume Boulant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Annexe : Règles de développement SALOME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cette annexe est un recueil de techniques de développement spécifiques
à l'environnement SALOME et utilisées pour la mise au point du module
XMED. Elles sont a priori utilisables pour d'autres contexte de
développement dans SALOME.

.. TODO: récupérer les fonctions génériques de VisuGUI_Tools.cxx
.. TODO: récupérer les fonctions génériques de SMESGGUI_utils.cxx

.. contents:: Sommaire
   :local:
   :backlinks: none

Récupérer la sélection dans l'arbre d'étude
===========================================

Dans une classe dérivée de ``SalomeApp_Module``, on peut utiliser un
code de la forme suivante:

.. code-block:: cpp
 
 #include <SALOME_ListIO.hxx>
 #include <LightApp_SelectionMgr.h>
 #include <SALOME_ListIteratorOfListIO.hxx>
 #include <SALOME_InteractiveObject.hxx>
 
 // ...

 // Get the selected object in the study (SObject)
 LightApp_SelectionMgr* aSelectionMgr = this->getApp()->selectionMgr();
 SALOME_ListIO aListIO;
 aSelectionMgr->selectedObjects(aListIO);

 // Analyse the selection. There can be more than one item.
 SALOME_ListIteratorOfListIO It (aListIO);
 for (; It.More(); It.Next()) {
   Handle(SALOME_InteractiveObject) anIO = It.Value();
   SALOMEDS::SObject_var aSObject = aStudy->FindObjectID(anIO->getEntry());

   // Check if the selected object is relevant for the operation
   // ...

   // Process the SObject if it's relevant
   // ...
   
 }

On peut noter qu'une variable ``aStudy`` est requise. Elle représente
l'étude SALOME sur laquelle s'oppère la sélection. L'étude active
(impliquée dans la sélection) peut être obtenue au moyen du
gestionnaire d'étude (voir :ref:`ci-dessous <salomedoc_getActiveStudy>`).

Réglage du curseur graphique
============================

Dans le cas où le traitement est long, il peut être intéressant
d'encadrer l'opération par un vérouillage du curseur de sélection:

.. code-block:: cpp

 QApplication::setOverrideCursor(Qt::WaitCursor);

 // Do the job
 // ...

 QApplication::restoreOverrideCursor();


Les variables pour la gestion de l'étude
========================================

Les variables CORBA
-------------------

Les variables CORBA comme le serveur de noms (naming service) et le
gestionnaire de cycle de vie des objets (LifeCycleCRORBA) sont
fréquement utilisés. Dans le contexte d'une application (classe de
type ``SalomeApp_Module``), il est possible de récupérer simplement
des instances de ces variables par les commandes suivantes:

.. code-block:: cpp

   #include <SalomeApp_Application.h>
   #include <SALOME_NamingService.hxx>

   CORBA::ORB_var               orb           = SalomeApp_Application::orb();
   SALOMEDSClient_StudyManager* studyMgr      = SalomeApp_Application::studyMgr();
   SALOME_NamingService*        namingService = SalomeApp_Application::namingService();
   SALOME_LifeCycleCORBA*       lcc           = SalomeApp_Application::lcc();

Pour un usage en dehors de l'application graphique (par exemple au
niveau du container), l'orb peut être obtenu par les mécanismes
standard d'omniORB:

.. code-block:: cpp

   CORBA::ORB_var orb = CORBA::ORB_init(0,0);

L'orb est par exemple utile à récupérer pour la sérialisation des
objets CORBA et la manipulation des références sous forme de chaîne de
caractères:

.. code-block:: cpp

   // We suppose here that we have a CORBA object reference (object of
   // type *_ptr or *_var), for example a SALOME_MED::MED object.
   SALOME_MED::MED_ptr medObj = ... // anything to get this object  
   QString medIOR = orb->object_to_string(medObj);

   SALOME_MED::MED_ptr anOtherRefToMedObj = orb->string_to_object(medIOR)

.. note: this serialization can be used to communicate between a GUI
   and a component in a container, or between the C++ context and the
   python context.

.. _salomedoc_getActiveStudy:

Récupérer l'étude active
------------------------

Le concept d'étude active est un concept GUI. Il désigne l'étude en
cours d'usage au niveau de l'interface graphique. 

.. note: Pour rappel, l'étude est un objet CORBA de type
   ``SALOMEDS::Study`` qui héberge physiquement les ``SObject``
   pointant vers les données. L'arbre d'étude ("Object browser") est
   une représentation graphique de cet objet.

L'étude active peut être obtenue au moyen du gestionnaire
d'étude. Dans le corps d'une classe de type ``SalomeApp_Module``, ceci
peut se faire par un code de la forme suivante:

.. code-block:: cpp

  #include "SALOMEconfig.h"
  #include CORBA_SERVER_HEADER(SALOMEDS)
  #include <SalomeApp_Application.h>
  #include <SALOME_NamingService.hxx>

  // ...

  // Get the study id of the active study
  SalomeApp_Study* appStudy = dynamic_cast<SalomeApp_Study*> (this->getApp()->activeStudy());
  _PTR(Study) aCStudy = appStudy->studyDS();
  int aStudyID = aCStudy->StudyId();

  // Then get the study manager
  SALOME_NamingService *aNamingService = SalomeApp_Application::namingService();
  CORBA::Object_var anObject = aNamingService->Resolve("/myStudyManager");
  SALOMEDS::StudyManager_var aStudyManager = SALOMEDS::StudyManager::_narrow(anObject);

  // Finally, request the study manager for the study (SALOMEDS::Study)
  SALOMEDS::Study_var aStudy = aStudyManager->GetStudyByID(aStudyID);


Communiquer avec la console python
==================================

La console python désigne l'interpréteur embarqué dans l'interface
graphique de SALOME (GUI). Elle est également désignée comme
l'interface textuelle de SALOME (TUI) car elle permet de piloter
SALOME au moyen de commandes en syntaxe python.

Le paragraphe montre comment communiquer avec cette interface texte
depuis le contexte C++ de l'interface graphique, en particulier pour
déclencher l'exécution de commandes.

Le code se situe donc au sein d'une classe de type
``SalomeApp_Module`` (de laquelle hérite la partie graphique d'un
module SALOME, de nom ``<MODULE_NAME>GUI``). Cette classe possède une
méthode ``getApp()`` par laquelle on peut récupérer une instance de la
console python embarquée (this->getApp()->pythonConsole()).

Le code suivant illustre l'envoie d'une commande python par ce
mécanisme. Dans cette exemple, on défini une variable ``id`` dans la
console python comme l'identifiant de l'étude active:

.. code-block:: cpp

   #include <PyConsole_Console.h>
   #include <QString>
   #include <QStringList>

   PyConsole_Console * pyConsole = getApp()->pythonConsole();

   QStringList commands;
   commands+="import salome";
   commands+="id=salome.myStudyId";
      
   QStringListIterator it(commands);
   while (it.hasNext()) {
       pyConsole->exec(it.next());
   }

Dans ce deuxième exemple, on cherche à reconstituer dans le contexte
de la console python un pointer vers un objet med instancié dans le
contexte C++ de l'application graphique. Pour cela, on communique la
référence de l'objet sous la forme sérialisé (IOR pour un objet
CORBA):

.. code-block:: cpp

   #include <PyConsole_Console.h>
   #include <QString>
   #include <QStringList>
   #include <SalomeApp_Application.h>

   // We suppose here that we have a CORBA object reference (object of
   // type *_ptr or *_var), for example a SALOME_MED::MED object.
   SALOME_MED::MED_ptr medObj = ... // anything to get this object  

   // Get the IOR of this object
   QString medIOR = SalomeApp_Application::orb()->object_to_string(medObj);

   PyConsole_Console * pyConsole = getApp()->pythonConsole();

   QStringList commands;
   commands+="import salome";
   commands+=QString("med=salome.orb.string_to_object(\"%1\")").arg(medIOR);
      
   QStringListIterator it(commands);
   while (it.hasNext()) {
       pyConsole->exec(it.next());
   }
