#!/usr/bin/env python

from MEDLoader import MEDLoader

import os
#filename = "madnex_field.med"
filename = "timeseries.med"
filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),filename)

# Read the source meshes
meshNames = MEDLoader.GetMeshNames(filepath)

# Set to True if the meshes and fields data must be loaded. Otherwise,
# only theire descriptions will be loaded.
READ_PHYSICAL_DATA=False

for meshName in meshNames:

    print "%s"%meshName

    # At this step, one can load the mesh of name meshName (but it is
    # not an obligation to continue to explore the metadata)
    meshDimRelToMax = 0 # 0 = no restriction
    if READ_PHYSICAL_DATA:
        mesh = MEDLoader.ReadUMeshFromFile(filepath,meshName,meshDimRelToMax)
    # Note that the read function required the parameter
    # meshDimRelToMax. This parameter discreminates the meshdim you
    # are interested to relatively to the maximal dimension of cells
    # contained in the mesh in file (then its value could be 0, -1, -2
    # or -3 depending on the max dimension of the mesh. 0 means "no
    # restriction".

    # Read the names of the fields that rely on this mesh
    fieldNames = MEDLoader.GetAllFieldNamesOnMesh(filepath,meshName)

    for fieldName in fieldNames:

        print "  %s"%fieldName
        
        # A field name could identify several MEDCoupling fields, that
        # differ by their spatial discretization on the mesh (values on
        # cells, values on nodes, ...). This spatial discretization is
        # specified by the TypeOfField that is an integer value in this
        # list:
        # 0 = ON_CELLS 	
        # 1 = ON_NODES 	
        # 2 = ON_GAUSS_PT 	
        # 3 = ON_GAUSS_NE
        #
        # As a consequence, before loading values of a field, we have
        # to determine the types of spatial discretization defined for
        # this field and to chooose one.

        listOfTypes = MEDLoader.GetTypesOfField(filepath,meshName,fieldName)
        for typeOfDiscretization in listOfTypes:
            print "    %s"%typeOfDiscretization

            # Then, we can get the iterations associated to this field on
            # this type of spatial discretization:
            fieldIterations = MEDLoader.GetFieldIterations(typeOfDiscretization,
                                                           filepath,
                                                           meshName,
                                                           fieldName)

            # Then, we can access to the physical data for each
            # iteration of this field
            for fieldIteration in fieldIterations:
                itNumber = fieldIteration[0]
                itOrder  = fieldIteration[1]
                print "      (%s,%s)"%(itNumber,itOrder)
                
                if READ_PHYSICAL_DATA:
                    medCouplingField = MEDLoader.ReadField(typeOfDiscretization,
                                                           filepath,
                                                           meshName,
                                                           meshDimRelToMax,
                                                           fieldName,
                                                           itNumber,
                                                           itOrder)
                    print medCouplingField
