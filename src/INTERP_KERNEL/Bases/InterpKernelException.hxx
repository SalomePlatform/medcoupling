// Copyright (C) 2007-2026  CEA, EDF
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

#ifndef __INTERPKERNELEXCEPTION_HXX__
#define __INTERPKERNELEXCEPTION_HXX__

#include "INTERPKERNELDefines.hxx"

#include <string>
#include <exception>
#include <sstream>

namespace INTERP_KERNEL
{
class Exception : public std::exception
{
   public:
    INTERPKERNEL_EXPORT Exception() noexcept;
    INTERPKERNEL_EXPORT Exception(const Exception &other) noexcept;
    INTERPKERNEL_EXPORT Exception(const char *reason);
    INTERPKERNEL_EXPORT Exception(const std::string &reason);
    INTERPKERNEL_EXPORT Exception(const char *reason, const char *file, int line);
    INTERPKERNEL_EXPORT ~Exception();
    INTERPKERNEL_EXPORT const char *what() const noexcept override;

   protected:
    std::string _reason;
};
}  // namespace INTERP_KERNEL

#define THROW_IK_EXCEPTION(text)                   \
    {                                              \
        std::ostringstream oss;                    \
        oss << text;                               \
        throw INTERP_KERNEL::Exception(oss.str()); \
    }

#endif
