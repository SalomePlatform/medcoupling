#include "RemapperTest.hxx"
#include "Remapper.hxx"

#include <iostream>
#include <vector>

namespace INTERP_TEST
{


  void RemapperTest::setUp() 
  {
  }

 
  void RemapperTest::tearDown() 
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
  void RemapperTest::test_Remapper() {
		string sourcename=getenv("MED_ROOT_DIR");
		sourcename +="/share/salome/resources/MedFiles/square1.med";
		MEDMEM::MESH source_mesh (MED_DRIVER,sourcename,"Mesh_2");

		string targetname=getenv("MED_ROOT_DIR");
		targetname +="/share/salome/resources/MedFiles/square2.med";
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


		INTERP_KERNEL::Remapper remapper;
		remapper.prepare(source_mesh,target_mesh);
		remapper.transfer(source_field,target_field);

		MEDMEM::FIELD<double> source_areas(&source_support,1);
		MEDMEM::FIELD<double> target_areas(&target_support,1);
		source_areas.getArea();
		target_areas.getArea();
		absField(source_areas); //absolute value
		absField(target_areas); //absolute value


		//target square is in reverse order as compared to initial square
		double source_integral=source_field.normL2(1,&source_areas);
		double target_integral=target_field.normL2(1,&target_areas);
		
		CPPUNIT_ASSERT_DOUBLES_EQUAL(source_integral,target_integral,1e-10);
		
  }

		void RemapperTest::absField(MEDMEM::FIELD<double>& field)
			{
				double* areas=const_cast<double*>(field.getValue());
				for (int i=0; i< field.getNumberOfValues();i++)
					{
						areas[i]=fabs(areas[i]);
					}
			}

}
