#ifndef _UNITTESTSRESULT_HXX_
#define _UNITTESTSRESULT_HXX_

#include <fstream>
#include <cstdlib>

namespace INTERP_KERNEL
{
  static inline std::string getResultFile() 
    {
      std::string s = "/tmp/";
      s += std::getenv("USER");
      s += "/UnitTestsResult";
      return s;
    }

  std::string UnitTestsResult = getResultFile();
}

#endif
