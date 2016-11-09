
MEDCoupling,  multiprocessing
----------------------------

Cet exercice fait la supposition de Numpy Scipy sont correctement maîtrisés, sinon voir :ref:`medcouplingnumpyptr`.
On va faire simuler un traitement un peu long (ici l'interpolation d'un maillage de 64000 cells avec un autre de 64000 cells).
On va faire le traitement d'abord en séquentiel puis en parallèle pour exploiter les coeurs de notre CPU.
Nous allons utiliser le module ``multiprocessing`` pour cela.

Début de l'implémentation
~~~~~~~~~~~~~~~~~~~~~~~~~

Pour commencer l'exercice importer le module Python ``MEDCoupling``, ``MEDCouplingRemapper``, ``numpy``, ``scipy``, ``multiprocessing``
et ``datetime`` pour chronométrer : ::

	import MEDCoupling as mc
	import MEDCouplingRemapper as mr
	import multiprocessing as mp
	from datetime import datetime
	from scipy.sparse import csr_matrix

Créer un maillage cartésien régulier 3D avec 40 cells en X, Y et Z
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Créons un maillage cartésien 3D m de pas régulier entre 0. et 1. en X, Y et Z : ::

	nbCells=40
	arr=mc.DataArrayDouble(nbCells+1) ; arr.iota() ; arr/=nbCells
	m=mc.MEDCouplingCMesh() ; m.setCoords(arr,arr,arr)

Traduisons m en non structuré pour que le calcul soit plus long : ::

	m=m.buildUnstructured()

Créer une copie m2 de m translatée de la moitié du pas en X, Y et Z ::

	m2=m.deepCopy()
	t=mc.DataArrayDouble(3)
	t[:]=1/(2*float(nbCells))
	m2.translate(t.getValues())

Calculer séquentiellement la matrice d'interpolation M de la projection entre m et m2 en P0P0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m sera considéré comme le maillage source et m2 sera considéré comme le maillage cible.
Profitons en pour chronométrer le temps necessaire pour le traitement séquentiel.
Utilisons ``MEDCouplingRemapper`` pour cela. ::

	remap=mr.MEDCouplingRemapper()
	strt=datetime.now()
	assert(remap.prepare(m,m2,"P0P0")==1)
	print "time in sequential : %s"%(str(datetime.now()-strt))

Stockons la sparse matrix scipy dans ``matSeq``. ::

	matSeq=remap.getCrudeCSRMatrix()

Calculer cette même matrice M en parallèle avec multiprocessing.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Commencons par récupérer le nombre de coeur de notre machine. ::

	nbProc=mp.cpu_count()

L'idée est de faire une méthode ``work`` prenant un tuple de longueur 2.
Le premier élément du tuple est la partie du maillage ``m2`` considérée. Le 2eme élément donne la correspondance entre les cells id de ``m2Part`` les cells id de ``m2``.

L'idée est d'interpoler ``m`` avec ``m2Part``.

On récupèrera ensuite la matrice sparse ``myMat`` issue de ``m`` avec ``m2Part``.
Ensuite l'idée et de générer une matrice sparse ``mat2`` à partir de ``myMat`` avec les ids globaux de ``m2``. ::

	def work(inp):
            m2Part,partToGlob=inp
	    myRemap=mr.MEDCouplingRemapper()
	    assert(myRemap.prepare(m,m2Part,"P0P0")==1)
	    myMat=myRemap.getCrudeCSRMatrix()
	    a=mc.DataArrayInt.Range(s.start,s.stop,s.step)
	    indptrnew=mc.DataArrayInt(m2.getNumberOfCells())
	    indptrnew.fillWithZero()
	    d=mc.DataArrayInt(myMat.indptr).deltaShiftIndex()
	    indptrnew[partToGlob]=d
	    indptrnew.computeOffsetsFull()
	    mat2=csr_matrix( (myMat.data,myMat.indices,indptrnew.toNumPyArray()), shape=(m2.getNumberOfCells(),m.getNumberOfCells()))
	    return mat2

Il s'agit désormais de faire la liste des inputs à donner aux ``nbProc`` instances de ``work`` qui seront exécutés en parallèle.
Appelons cette liste python ``workToDo`` qui sera de longueur ``nbProc``.
On peut se faire aider de ``mc.DataArray.GetSlice``. ::

        workToDo=[]
        for i in xrange(nbProc):
              s=mc.DataArray.GetSlice(slice(0,m2.getNumberOfCells(),1),i,nbProc)
              part=m2[s]
              partToGlob=mc.DataArrayInt.Range(s.start,s.stop,s.step)
              workToDo.append((part,partToGlob))
              pass

On est maintenant prêt pour faire travailler chacun des coeurs indépendamment. Pour ce faire, on crée un ``mp.Pool`` et on assigne à chaque worker le travail ``work`` avec autant de worker que de coeurs. Et chronométrons le tout ! ::


	strt=datetime.now()
	pool = mp.Pool()
	asyncResult = pool.map_async(work,workToDo)
	subMatrices = asyncResult.get()
	print "time in parallel (x%d) : %s"%(nbProc,str(datetime.now()-strt))

.. note:: A noter la magie ! On a transféré entre le process maitre et les process esclave sans même s'en rendre compte les maillages et les DataArrayInt contenus dans ``workToDo`` !
	  Merci à la pickelisation des objets MEDCoupling :)

Vérfication
~~~~~~~~~~~

Vérifions que les matrices sont les mêmes ! Sommons ``subMatrices`` (``matPar``) et regardons le nombre de non zéros de la différence entre la ``matPar`` et ``matSeq``. ::

	matPar = sum(subMatrices)
	matDelta=matSeq-matPar
	assert(matDelta.nnz==0)
