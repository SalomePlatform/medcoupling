#ifndef RENUMBERINGFACTORY_HXX_
#define RENUMBERINGFACTORY_HXX_

#include <string>
#include "RENUMBER_Renumbering.hxx"

namespace MED_RENUMBER
{
  Renumbering* RenumberingFactory(const std::string& s);
}

#endif /*RENUMBERINGFACTORY_HXX_*/
