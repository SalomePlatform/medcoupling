#include "RenumberingFactory.hxx"
#include "RENUMBER_Renumbering.hxx"
#ifdef ENABLE_METIS
#include "RENUMBER_METISRenumbering.hxx"
#endif
#ifdef ENABLE_BOOST
#include "RENUMBER_BOOSTRenumbering.hxx"
#endif

#include <iostream>

using namespace std;

namespace MED_RENUMBER
{
  Renumbering* RenumberingFactory(const string &s)
  {
#ifdef ENABLE_METIS
#ifdef ENABLE_BOOST
    if (s=="METIS")
      {
        return new METISRenumbering;
      }
    else if(s=="BOOST")
      {
        return new BOOSTRenumbering;
      }
    else 
      {
        std::cerr << "The method has to be METIS or BOOST" << std::endl;
        return 0;
      }
#endif
#ifndef ENABLE_BOOST
    if (s=="METIS")
      {
        return new METISRenumbering;
      }
    else
      {
        std::cerr << "The method has to be METIS!" << std::endl;
        return 0;
      }
#endif
#endif
#ifndef ENABLE_METIS
#ifdef ENABLE_BOOST
    if (s=="BOOST")
      {
        return new BOOSTRenumbering;
      }
    else
      {
        std::cerr << "The method has to be BOOST!" << std::endl;
        return 0;
      }
#endif
#ifndef ENABLE_BOOST
    std::cerr << "Error, no method compiled" << std::endl;
    return 0;
#endif
#endif
  }
}
