// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDFileMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "InterpolationOptionsTest.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "Interpolation2D.txx"
#include "TestInterpKernelUtils.hxx"

#include <iostream>
#include <vector>

using namespace MEDCoupling;

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
    std::string sourcename=INTERP_TEST::getResourceFile("square1.med");
    MEDFileUMesh *source_mesh=MEDFileUMesh::New(sourcename.c_str(),"Mesh_2");

    std::string targetname=INTERP_TEST::getResourceFile("square2.med");
    MEDFileUMesh *target_mesh=MEDFileUMesh::New(targetname.c_str(),"Mesh_3");

    MEDCouplingUMesh *source_mesh_mc=source_mesh->getMeshAtLevel(0);
    MEDCouplingFieldDouble *source_field=MEDCouplingFieldDouble::New(ON_CELLS);
    source_field->setMesh(source_mesh_mc); source_mesh_mc->decrRef();
    DataArrayDouble *arr=DataArrayDouble::New(); arr->alloc(source_mesh_mc->getNumberOfCells(),1);
    source_field->setArray(arr); arr->decrRef();
    double *value=arr->getPointer();
    for(int i=0; i<source_mesh_mc->getNumberOfCells(); i++)
      value[i]=1.0;
    MEDCouplingUMesh *target_mesh_mc=target_mesh->getMeshAtLevel(0);
    MEDCouplingFieldDouble *target_field=MEDCouplingFieldDouble::New(ON_CELLS);
    target_field->setMesh(target_mesh_mc); target_mesh_mc->decrRef();
    arr=DataArrayDouble::New(); arr->alloc(target_mesh_mc->getNumberOfCells(),1);
    target_field->setArray(arr); arr->decrRef();
    double* targetvalue=arr->getPointer();
    for(int i=0; i<target_mesh_mc->getNumberOfCells(); i++)
      targetvalue[i]=0.0;
    // Ok at this point we have our mesh in MED-Memory format.
    // Go to wrap med_source_mesh and med_target_mesh.
    MEDCouplingNormalizedUnstructuredMesh<2,2> wrap_source_mesh(source_mesh_mc);
    MEDCouplingNormalizedUnstructuredMesh<2,2> wrap_target_mesh(target_mesh_mc);
    // Go for interpolation...
    INTERP_KERNEL::Interpolation2D myInterpolator;
    //optionnal call to parametrize your interpolation. First precision, tracelevel, intersector wanted.
    myInterpolator.setPrecision(1e-7);
    myInterpolator.setPrintLevel(1);
    source_mesh->decrRef();
    source_field->decrRef();
    target_field->decrRef();
    target_mesh->decrRef();
  }
}
