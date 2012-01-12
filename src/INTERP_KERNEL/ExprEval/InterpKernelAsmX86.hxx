// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __INTERPKERNELASMX86_HXX__
#define __INTERPKERNELASMX86_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"

#include <vector>
#include <string>

namespace INTERP_KERNEL
{
  class AsmX86
  {
  public:
    std::vector<char> convertIntoMachineLangage(const std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    char *copyToExecMemZone(const std::vector<char>& ml, unsigned& offset) const;
  private:
    void convertOneInstructionInML(const std::string& inst, std::vector<char>& ml) const throw(INTERP_KERNEL::Exception);
  private:
    static void convertMov(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertPush(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertPop(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFld(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFaddp(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFsubp(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFmulp(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFdivp(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFcos(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFsin(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFabs(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFchs(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFsqrt(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertSub(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertAdd(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertRet(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertLeave(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertMovsd(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void convertFst(const std::string& inst, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    //
    static void convertMovToEsp(const std::string& inst1, const std::string& inst2, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
    static void appendAddress(const std::string& addr, int nbOfByte, std::vector<char>& ml) throw(INTERP_KERNEL::Exception);
  private:
    static const int NB_OF_OPS=19;
    static const char *OPS[NB_OF_OPS];
  };
}

#endif
