// Copyright (C) 2021-2024  CEA, EDF
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

#include "MEDMESHConverterUtilities.hxx"

#include "libmesh5.hxx"

namespace MeshFormat
{
    bool isMeshExtensionCorrect( const std::string& fileName )
    {
        std::string ext;
        {
            auto pos(fileName.find_last_of('.'));
            if(pos != std::string::npos)
                ext = fileName.substr(pos);
        }
        switch ( ext.size() ) {
        case 5:
            return ( ext == ".mesh" || ext == ".solb" );
        case 6:
            return ( ext == ".meshb" );
        case 4:
            return ( ext == ".sol" );
        }
        return false;
    }


    Localizer::Localizer()
    {
        _locale = setlocale(LC_NUMERIC, NULL);
        setlocale(LC_NUMERIC, "C");
    }

    Localizer::~Localizer()
    {
        setlocale(LC_NUMERIC, _locale.c_str() );
    }
}


