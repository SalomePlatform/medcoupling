//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef MEDSPLITTER_MESHCOLLECTIONMEDASCIIDRIVER_H
#define MEDSPLITTER_MESHCOLLECTIONMEDASCIIDRIVER_H


#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

/*!reads a distributed field
 * 
 * \param fields vector of fields (one field per subdomain)
 * \param fieldname name of the field
 * \param itnumber number of iteration
 * \param ordernumber internal number inside the iteration
 * */
template <class T>
void MESHCollectionMedAsciiDriver::_readFields(vector<MEDMEM::FIELD<T>* >& fields,char* fieldname, int itnumber, int ordernumber)
{
  for (int i=0; i<_collection->getMesh().size(); i++)
  {
    char filename[256];
    strcpy(filename,_filename[i].c_str());
    cout << "maillage : " << filename << " champ : " << fieldname << endl;
    //    MEDMEM::FIELD<T>* field = new MEDMEM::FIELD<T>(MEDMEM::MED_DRIVER,filename,fieldname,itnumber,ordernumber);
    fields.push_back (new MEDMEM::FIELD<T>(MEDMEM::MED_DRIVER,filename,fieldname,itnumber,ordernumber));
  }
}


/*!writes a distributed field
 * 
 * \param fields vector of fields (one field per subdomain)
 * \param fieldname name of the field
 * */
template <class T>
void MESHCollectionMedAsciiDriver::_writeFields(vector<MEDMEM::FIELD<T>* >& fields,char* fieldname)
{
  for (int i=0; i<_collection->getMesh().size(); i++)
  {
    char filename[256];
    strcpy(filename,_filename[i].c_str());
    int driverid = fields[i]->addDriver(MEDMEM::MED_DRIVER, filename, fieldname);
    fields[i]->write(driverid);    
  }
}

#endif
