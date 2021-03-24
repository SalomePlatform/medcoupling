// Copyright (C) 2021  CEA/DEN, EDF R&D
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

#ifndef MEDMESHCONVERTERUTILITIES_HXX
#define MEDMESHCONVERTERUTILITIES_HXX

#include "MEDLoaderDefines.hxx"

# include <string>
# include <sstream>
#include <iostream>
# include <vector>
#include <cstring>



#define THROW_IK_EXCEPTION(text)                        \
{                                                     \
    std::ostringstream oss; oss << text;                \
    throw INTERP_KERNEL::Exception(oss.str().c_str());  \
}

namespace MeshFormat
{
bool isMeshExtensionCorrect( const std::string& fileName );

//~void closeMeshFile(const int gmfMeshID);

enum Status {
    DRS_OK,
    DRS_EMPTY,          // a file contains no mesh with the given name
    DRS_WARN_RENUMBER,  // a file has overlapped ranges of element numbers,
    // so the numbers from the file are ignored
    DRS_WARN_SKIP_ELEM, // some elements were skipped due to incorrect file data
    DRS_WARN_DESCENDING, // some elements were skipped due to descending connectivity
    DRS_FAIL            // general failure (exception etc.)
};
/*!
 * \brief Class to generate string from any type
 */
class Comment : public std::string
{
    std::ostringstream _s ;

public :

    Comment():std::string("") {}

    Comment(const Comment& c):std::string() {
        _s << c.c_str() ;
        this->std::string::operator=( _s.str() );
    }

    Comment & operator=(const Comment& c) {
        _s << c.c_str() ;
        this->std::string::operator=( _s.str() );
        return *this;
    }

    template <class T>
    Comment( const T &anything ) {
        _s << anything ;
        this->std::string::operator=( _s.str() );
    }

    template <class T>
    Comment & operator<<( const T &anything ) {
        _s << anything ;
        this->std::string::operator=( _s.str() );
        return *this ;
    }

    operator char*() const {
        return (char*)c_str();
    }

    std::ostream& Stream() {
        return _s;
    }
};
/*!
   * \brief Enforce freeing memory allocated by std::vector
   */
template <class TVECTOR>
void FreeVector(TVECTOR& vec)
{
    TVECTOR v2;
    vec.swap( v2 );
}
template <class TVECTOR>
void CompactVector(TVECTOR& vec)
{
    TVECTOR v2( vec );
    vec.swap( v2 );
}

class Localizer
{
    std::string _locale;
public:
	MEDLOADER_EXPORT Localizer();
	MEDLOADER_EXPORT ~Localizer();
};
}


#endif // MEDMESHCONVERTERUTILITIES_HXX
