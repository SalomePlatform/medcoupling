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

#include "ParaMEDMEMTest.hxx"
#include "TestInterpKernelUtils.hxx"
#include <cppunit/TestAssert.h>

#include <sstream>
#include <cmath>
#include <list>
#include <stdexcept>
#include <stdlib.h>

#ifndef WIN32
#include <unistd.h>
#endif



//================================================================================
/*!
 * \brief Returns writable temporary directory
 * \return full path to the temporary directory
 */
//================================================================================

std::string ParaMEDMEMTest::getTmpDirectory()
{
  std::string path;

  std::list<std::string> dirs;
  if ( getenv("TMP") )    dirs.push_back( getenv("TMP" ));
  if ( getenv("TMPDIR") ) dirs.push_back( getenv("TMPDIR" ));
  dirs.push_back( "/tmp" );

  std::string tmpd = "";
  for ( std::list<std::string>::iterator dir = dirs.begin(); dir != dirs.end() && tmpd == "" ; ++dir ) {
    if ( access( dir->data(), W_OK ) == 0 ) {
      tmpd = dir->data();
    }
  }

  if ( tmpd == "" )
    throw std::runtime_error("Can't find writable temporary directory. Set TMP environment variable");

  return tmpd;
}

//================================================================================
/*!
 * \brief Creates a copy of source file (if source file is specified)
 * in the temporary directory and returns a path to the tmp file
 *
 * \param tmpfile name of the temporary file (without path)
 * \param srcfile source file
 * \return path to the temporary file
 */
//================================================================================
std::string ParaMEDMEMTest::makeTmpFile( const std::string& tmpfile, const std::string& srcfile )
{
  std::string tmpf = getTmpDirectory() + "/" + tmpfile;
  if ( srcfile != "" ) {
    std::string cmd  = "cp " + srcfile + " " + tmpf + " ; chmod +w " + tmpf;
    system( cmd.c_str() );
  }
  return tmpf;
}


/*!
 *  Tool to remove temporary files.
 *  Allows automatique removal of temporary files in case of test failure.
 */
ParaMEDMEMTest_TmpFilesRemover::~ParaMEDMEMTest_TmpFilesRemover()
{
  std::set<std::string>::iterator it = myTmpFiles.begin();
  for (; it != myTmpFiles.end(); it++) {
    if (access((*it).data(), F_OK) == 0)
      remove((*it).data());
  }
  myTmpFiles.clear();
  //cout << "~ParaMEDMEMTest_TmpFilesRemover()" << endl;
}

bool ParaMEDMEMTest_TmpFilesRemover::Register(const std::string theTmpFile)
{
  return (myTmpFiles.insert(theTmpFile)).second;
}
