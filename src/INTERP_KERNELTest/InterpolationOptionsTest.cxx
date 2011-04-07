//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "InterpolationOptionsTest.hxx"
#include "MEDNormalizedUnstructuredMesh.txx"
#include "Interpolation2D.txx"
#include "TestInterpKernelUtils.hxx"
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
  void InterpolationOptionsTest::test_InterpolationOptions() 
  {
    string sourcename=INTERP_TEST::getResourceFile("square1.med");
    MEDMEM::MESH *source_mesh=new MEDMEM::MESH(MED_DRIVER,sourcename,"Mesh_2");

    string targetname=INTERP_TEST::getResourceFile("square2.med");
    MEDMEM::MESH *target_mesh=new MEDMEM::MESH(MED_DRIVER,targetname,"Mesh_3");

    const MEDMEM::SUPPORT *source_support=source_mesh->getSupportOnAll(MED_EN::MED_CELL);
    MEDMEM::FIELD<double> *source_field=new FIELD<double>(source_support,1);
    double* value=const_cast<double*>(source_field->getValue());
    for (int i=0; i<source_support->getNumberOfElements(MED_EN::MEDMEM_ALL_ELEMENTS); i++)
      value[i]=1.0;
    const MEDMEM::SUPPORT *target_support=target_mesh->getSupportOnAll(MED_EN::MED_CELL);
    MEDMEM::FIELD<double> *target_field=new FIELD<double>(target_support,1);
    double* targetvalue=const_cast<double*>(target_field->getValue());
    for (int i=0; i<target_support->getNumberOfElements(MED_EN::MEDMEM_ALL_ELEMENTS); i++)
      targetvalue[i]=0.0;
    // Ok at this point we have our mesh in MED-Memory format.
    // Go to wrap med_source_mesh and med_target_mesh.
    MEDNormalizedUnstructuredMesh<2,2> wrap_source_mesh(source_mesh);
    MEDNormalizedUnstructuredMesh<2,2> wrap_target_mesh(target_mesh);
    // Go for interpolation...
    INTERP_KERNEL::Interpolation2D myInterpolator;
    //optionnal call to parametrize your interpolation. First precision, tracelevel, intersector wanted.
    myInterpolator.setPrecision(1e-7);
    myInterpolator.setPrintLevel(1);
    source_mesh->removeReference();
    source_field->removeReference();
    target_field->removeReference();
    target_mesh->removeReference();
  }
}
