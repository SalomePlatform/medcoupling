#include "InterpolationOptionsTest.hxx"
#include "MEDNormalizedUnstructuredMesh.hxx"
#include "Interpolation2D.hxx"
#include <iostream>
#include <vector>

namespace INTERP_TEST
{


  void InterpolationOptionsTest::setUp() 
  {
  }

 
  void InterpolationOptionsTest::tearDown() 
  {
  }

  /**
   * Test that creates a tree in 2D and check that 
   * the results are correct in three
   * cases :
   * a non matching search
   * a standard case
   * a bbox overlapping the bboxes of the tree
   */
  void InterpolationOptionsTest::test_InterpolationOptions() {
		string sourcename=getenv("MED_ROOT_DIR");
		sourcename +="/share/salome/resources/med/square1.med";
		MEDMEM::MESH source_mesh (MED_DRIVER,sourcename,"Mesh_2");

		string targetname=getenv("MED_ROOT_DIR");
		targetname +="/share/salome/resources/med/square2.med";
		MEDMEM::MESH target_mesh (MED_DRIVER,targetname,"Mesh_3");

		MEDMEM::SUPPORT source_support(&source_mesh,"on All support");
		MEDMEM::FIELD<double> source_field(&source_support,1);
		double* value=const_cast<double*>(source_field.getValue());
		for (int i=0; i<source_support.getNumberOfElements(MED_EN::MED_ALL_ELEMENTS); i++)
			value[i]=1.0;
		
		MEDMEM::SUPPORT target_support(&target_mesh,"on All support");
		MEDMEM::FIELD<double> target_field(&target_support,1);
		double* targetvalue=const_cast<double*>(target_field.getValue());
		for (int i=0; i<target_support.getNumberOfElements(MED_EN::MED_ALL_ELEMENTS); i++)
			targetvalue[i]=0.0;

		// Ok at this point we have our mesh in MED-Memory format.
		// Go to wrap med_source_mesh and med_target_mesh.
		MEDNormalizedUnstructuredMesh<2,2> wrap_source_mesh(&source_mesh);
		MEDNormalizedUnstructuredMesh<2,2> wrap_target_mesh(&target_mesh);
		// Go for interpolation...
		INTERP_KERNEL::Interpolation2D myInterpolator;
		//optionnal call to parametrize your interpolation. First precision, tracelevel, intersector wanted.
		myInterpolator.setPrecision(1e-7);
		myInterpolator.setPrintLevel(1);
  }


}
