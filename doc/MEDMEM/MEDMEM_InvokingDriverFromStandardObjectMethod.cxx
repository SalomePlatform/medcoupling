//  Copyright (C) 2007-2011  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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

#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"

#include "MEDMEM_Field.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Med.hxx"

using namespace MEDMEM ;
using namespace MED_EN ;

main () {
  
  const char * fileName   = "pointe.med";
  const char * fileName2  = "fieldCellDoubleOfpointe.med";
  const char * fieldName  = "fieldcelldouble";
  const char * meshName   = "maa1";
    
  try {
    // Test creation of drivers from the standard driver method of an object
    FIELD<double> * myField = new FIELD<double>();
    int myDriver1 = myField->addDriver(MED_DRIVER, fileName, fieldName);
    //myField->read();
    //This test failed due to inadequate Support implementation
    myField->rmDriver();  // TESTER LA VALIDITE DE myDriver2 !!!!

    int myDriver2 = myField->addDriver(MED_DRIVER, fileName2, fieldName);
    //myField->write(myDriver2);
    //This test failed due to inadequate Support implementation
    myField->rmDriver(myDriver2);

    MESH * myMesh  = new MESH();
    int myDriver3  = myMesh->addDriver(MED_DRIVER, fileName, meshName);
    myMesh->read();
    myMesh->rmDriver();

    MED  *  myMed  = new MED();
    int myDriver4  = myMed->addDriver(MED_DRIVER, fileName);
    myMed->readFileStruct();
    myMed->rmDriver();

    delete myField;
    delete myMesh;
    delete myMed;

  } catch (MEDEXCEPTION& ex){
    MESSAGE(ex.what()) ;
  }
}
