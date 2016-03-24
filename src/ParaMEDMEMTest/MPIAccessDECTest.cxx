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

#include "MPIAccessDECTest.hxx"
#include <cppunit/TestAssert.h>

#include <sstream>
#include <cmath>

#ifndef WIN32
#include <unistd.h>
#endif

using namespace std;



/*!
 *  Tool to remove temporary files.
 *  Allows automatique removal of temporary files in case of test failure.
 */
MPIAccessDECTest_TmpFilesRemover::~MPIAccessDECTest_TmpFilesRemover()
{
  set<string>::iterator it = myTmpFiles.begin();
  for (; it != myTmpFiles.end(); it++) {
    if (access((*it).data(), F_OK) == 0)
      remove((*it).data());
  }
  myTmpFiles.clear();
  //cout << "~MPIAccessTest_TmpFilesRemover()" << endl;
}

bool MPIAccessDECTest_TmpFilesRemover::Register(const string theTmpFile)
{
  return (myTmpFiles.insert(theTmpFile)).second;
}
