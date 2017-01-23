// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#include "MEDFileBlowStrEltUp.hxx"

using namespace MEDCoupling;

MEDFileBlowStrEltUp::MEDFileBlowStrEltUp(const MEDFileFields *fsOnlyOnSE, const MEDFileMeshes *ms, const MEDFileStructureElements *ses)
{
  if(!fsOnlyOnSE || !ms || !ses)
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp constructor : NULL input pointer !");
  _ms.takeRef(ms); _ses.takeRef(ses);
  std::vector< std::pair<std::string,std::string> > ps;
  fsOnlyOnSE->getMeshSENames(ps);
  std::size_t sz(ps.size());
  _elts.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const std::pair<std::string,std::string>& p(ps[i]);
      MCAuto<MEDFileFields> f(fsOnlyOnSE->partOfThisLyingOnSpecifiedMeshSEName(p.first,p.second));
      _elts[i]=f;
    }
}
