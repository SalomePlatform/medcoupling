.. Code Integration and Code Coupling documentation master file, created by sphinx-quickstart on Tue Apr 28 14:31:38 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

medcoupling user's manual
=========================

medcoupling gathers several powerful functionalities around the input and output data of field-physics-oriented simulation codes. Data manipulated by medcoupling are objects relative to fields for simulation codes.

medcoupling functionalities are accessible through python modules.

medcoupling functionalities can be split into 4 categories:
  #. :doc:`Data movement <data_movement>`: read/write from/to file, reduce, extract, duplicate, aggregate, compare, exchange data memory to memory across process (image of multifile to file)
  #. :doc:`Data analysis <data_analysis>`: extract/gather information in data to transform it in a condensate workable data (not necessarely field linked) for further use.
  #. :doc:`Data conversion <data_conversion>`: interpolate, project, repare, decimate, format convertion to make discuss simulation codes each other
  #. :doc:`Data optimization <data_optimization>` for simulation code : renumbering, partition for multiprocessor codes

First, this documentation introduces :doc:`fundamental concepts/objects <basic_concepts>`  of medcoupling for a better understanding of examples.

.. toctree::
   :maxdepth: 2
   :numbered:
   :hidden:

   basic_concepts
   data_movement
   data_analysis
   data_conversion
   data_optimization

