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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDLOADERBASE_HXX__
#define __MEDLOADERBASE_HXX__

#include "MEDLoaderDefines.hxx"
#include "InterpKernelException.hxx"

#include <string>

class MEDLOADER_EXPORT MEDLoaderBase
{
public:
  static int getStatusOfFile(const std::string& fileName);
  static char *buildEmptyString(int lgth);
  static void getDirAndBaseName(const std::string& fullName, std::string& dirName, std::string& baseName);
  static std::string getPathSep();
  static std::string joinPath(const std::string& dirName, const std::string& baseName);
  static std::string buildUnionUnit(const char *name, int nameLgth, const char *unit, int unitLgth);
  static void splitIntoNameAndUnit(const std::string& s, std::string& name, std::string& unit);
  static void strip(std::string& s);
  static void safeStrCpy(const char *src, int maxLgth, char *dest, int behaviour);
  static void safeStrCpy2(const char *src, int maxLgth, char *dest, int behaviour);
  static std::string buildStringFromFortran(const char *expr, int lgth);
  static void zipEqualConsChar(std::string& s, int minConsSmChar);
  static std::string zipString(const std::string& src, int sizeToRespect);
public:
  static const int EXIST_RW=0;
  static const int NOT_EXIST=1;
  static const int EXIST_RDONLY=2;
  static const int EXIST_WRONLY=3;
  static const int DIR_LOCKED=4;
  static const char WHITE_SPACES[];
};

#endif
