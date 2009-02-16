//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
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


    INTERP_KERNEL::Remapper remapper;
    remapper.prepare(source_mesh,target_mesh,"P0P0");
    remapper.transfer(source_field,target_field);

    MEDMEM::FIELD<double> *source_areas=source_mesh.getArea(&source_support);
    MEDMEM::FIELD<double> *target_areas=target_mesh.getArea(&target_support);
    absField(*source_areas); //absolute value
    absField(*target_areas); //absolute value

    //target square is in reverse order as compared to initial square
    double source_integral=source_field.normL2(1,source_areas);
    double target_integral=target_field.normL2(1,target_areas);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(source_integral,target_integral,1e-10);
    delete source_areas;
    delete target_areas;
    
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
