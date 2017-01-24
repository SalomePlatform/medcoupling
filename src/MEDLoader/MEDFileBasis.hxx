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

#ifndef __MEDFILEBASIS_HXX__
#define __MEDFILEBASIS_HXX__

#include "InterpKernelException.hxx"

#include <string>
#include <vector>

namespace MEDCoupling
{
  class MEDFileString
  {
  public:
    MEDFileString(int maxLgth);
    ~MEDFileString();
    void clear();
    void set(const char *s);
    char *getPointer() { return _content; }
    const char *getReprForWrite() const { return _content; }
    std::string getRepr() const;
  private:
    int _max_lgth;
    char *_content;
  };


  class MEDFileMultiString
  {
  public:
    MEDFileMultiString(int nbOfCompo, int maxLgthPerCompo);
    ~MEDFileMultiString();
    void set(int compoId, const char *s);
    const char *getReprForWrite() const;
    std::vector<std::string> getRepr() const;
    std::string getReprPerComp(int compId) const;
  private:
    int _nb_of_comp;
    int _max_lgth_per_comp;
    char *_content;
  };
}

#endif
