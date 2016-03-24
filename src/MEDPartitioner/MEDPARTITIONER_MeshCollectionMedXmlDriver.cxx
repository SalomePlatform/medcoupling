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

#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_MeshCollectionDriver.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_MeshCollectionMedXmlDriver.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"

#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#ifdef WIN32
#include <time.h>
#include <windows.h>
#endif


#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

using namespace MEDPARTITIONER;

/*!\class MeshCollectionMedXmlDriver
 *
 *\brief Driver for MED 3.2 files having Xml master files
 *
 * Driver for reading and writing distributed files
 * for which the master file is written in an Xml format compliant with
 * the MED 3.2 specification.
 * The reading and writing of the meshes and fields are apart :
 * the meshes must always be written/read before the fields. Reading/Writing fields
 * is optional and is done field after field. API for reading/writing fields
 * is written with a template so that FIELD<int> and FIELD<double>
 * can be conveniently handled.
 */

MeshCollectionMedXmlDriver::MeshCollectionMedXmlDriver(MeshCollection* collection):MeshCollectionDriver(collection)
{
}

/*!reads a MED File Xml Master File v>=2.3
 * and mounts the corresponding meshes in memory
 * the connect zones are created from the joints
 * 
 *\param filename Xml file containing the list of MED v2.3 files
 * */

int MeshCollectionMedXmlDriver::read(const char* filename, ParaDomainSelector* domainSelector)
{
  //distributed meshes
  int nbdomain;
  _master_filename=filename;

  //reading ascii master file
  try
    {
      //Setting up the Xml tree corresponding to filename
      xmlDocPtr master_doc=xmlParseFile(filename);

      if (!master_doc)
        throw INTERP_KERNEL::Exception("Xml Master File does not exist or is not compliant with Xml scheme");

      //number of domains
      xmlXPathContextPtr xpathCtx = xmlXPathNewContext(master_doc);
      xmlXPathObjectPtr xpathObj = xmlXPathEvalExpression(BAD_CAST "//splitting/subdomain", xpathCtx);
      if (xpathObj==0 || xpathObj->nodesetval->nodeNr ==0)
        throw INTERP_KERNEL::Exception("Xml Master File does not contain /MED/splitting/subdomain node");

      //as subdomain has only one property which is "number"
      //it suffices to take the content of its first child 
      const char* mystring = (const char*)xpathObj->nodesetval->nodeTab[0]->properties->children->content;
      sscanf(mystring, "%d", &nbdomain);

      //mesh name
      xmlXPathFreeObject(xpathObj);
      xpathObj = xmlXPathEvalExpression(BAD_CAST "//content/mesh", xpathCtx);
      if (xpathObj==0 || xpathObj->nodesetval->nodeNr ==0)
        throw INTERP_KERNEL::Exception("Xml Master File does not contain /MED/content/mesh node");
      _collection->setName( (const char*)xpathObj->nodesetval->nodeTab[0]->properties->children->content);

      //cout << "nb domain " << nbdomain << endl;
      MyGlobals::_File_Names.resize(nbdomain);
      MyGlobals::_Mesh_Names.resize(nbdomain);
      (_collection->getMesh()).resize(nbdomain);
      (_collection->getFaceMesh()).resize(nbdomain);
      (_collection->getCellFamilyIds()).resize(nbdomain);
      (_collection->getFaceFamilyIds()).resize(nbdomain);

      //retrieving the node which contains the file names
      const char filechar[]="//files/subfile";
      xmlXPathFreeObject(xpathObj);
      xpathObj = xmlXPathEvalExpression(BAD_CAST filechar, xpathCtx);
      if (xpathObj==0 || xpathObj->nodesetval->nodeNr ==0)
        throw INTERP_KERNEL::Exception("Xml Master File does not contain /MED/files/subfile nodes");
      int nbfiles = xpathObj->nodesetval ->nodeNr;
    
      for (int i=0; i<nbfiles;i++)
        {
          //reading information about the domain
          std::string host;
          //reading file names 
          std::ostringstream name_search_string;
          name_search_string<<"//files/subfile[@id=\""<<i+1<<"\"]/name";
          xmlXPathObjectPtr xpathObjfilename =
            xmlXPathEvalExpression(BAD_CAST name_search_string.str().c_str(),xpathCtx);
          if (xpathObjfilename->nodesetval ==0)
            throw INTERP_KERNEL::Exception("Error retrieving a file name from subfile of Xml Master File");
          MyGlobals::_File_Names[i]=(const char*)xpathObjfilename->nodesetval->nodeTab[0]->children->content;

          //reading the local mesh names
          std::ostringstream mesh_search_string;
          mesh_search_string<<"//mapping/mesh/chunk[@subdomain=\""<<i+1<<"\"]/name";

          xmlXPathObjectPtr xpathMeshObj = xmlXPathEvalExpression(BAD_CAST mesh_search_string.str().c_str(),xpathCtx);
          if (xpathMeshObj->nodesetval ==0)
            throw INTERP_KERNEL::Exception("Error retrieving mesh name from chunk of Xml Master File");
          MyGlobals::_Mesh_Names[i]=(const char*)xpathMeshObj->nodesetval->nodeTab[0]->children->content;

          if ( !domainSelector || domainSelector->isMyDomain(i))
            readSubdomain(i);
          xmlXPathFreeObject(xpathObjfilename);
          xmlXPathFreeObject(xpathMeshObj);
        } //loop on domains

      //LIBXML cleanup
      xmlXPathFreeObject(xpathObj); 
      xmlXPathFreeContext(xpathCtx); 
      xmlFreeDoc(master_doc); 

    } //of try
  catch(...)
    {
      throw INTERP_KERNEL::Exception("I/O error reading parallel MED file");
    }

  ParallelTopology* aPT = new ParallelTopology(_collection->getMesh());
  //creation of topology from mesh and connect zones
  if ( _collection->isParallelMode() )
    {
      //to know nb of cells on each proc to compute global cell ids from locally global
      domainSelector->gatherNbOf(_collection->getMesh());
    }
  _collection->setTopology(aPT, true);
  _collection->setDomainNames(_collection->getName());
  return 0;
}


/*! writes the collection of meshes in a 
 * MED v2.3 Xml file
 * with the connect zones being written as joints
 * \param filename name of the Xml file containing the meshes description
 */
void MeshCollectionMedXmlDriver::write(const char* filename, ParaDomainSelector* domainSelector) const
{
  xmlDocPtr master_doc = 0;
  xmlNodePtr root_node = 0, node, node2;
  char buff[256];

  //Creating the Xml document
  master_doc = xmlNewDoc(BAD_CAST "1.0");
  root_node = xmlNewNode(0, BAD_CAST "root");
  xmlDocSetRootElement(master_doc,root_node);

  //Creating child nodes
  // Version tag
  node = xmlNewChild(root_node, 0, BAD_CAST "version",0);
  xmlNewProp(node, BAD_CAST "maj", BAD_CAST "2");
  xmlNewProp(node, BAD_CAST "min", BAD_CAST "3");
  xmlNewProp(node, BAD_CAST "ver", BAD_CAST "1");

  //Description tag
  time_t present; 
  char date[20];
#ifndef WIN32
  time( &present);
  struct tm *time_asc = localtime(&present);
  sprintf(date,"%02d%02d%02d",time_asc->tm_year
          ,time_asc->tm_mon+1
          ,time_asc->tm_mday);
#else
  SYSTEMTIME    st;
  GetLocalTime ( &st );
  sprintf(date,"%02d%02d%02d",
          st.wYear
          ,st.wMonth
          ,st.wDay);
#endif

  node = xmlNewChild(root_node,0, BAD_CAST "description",0);

  xmlNewProp(node, BAD_CAST "what", BAD_CAST _collection->getDescription().c_str());
  xmlNewProp(node, BAD_CAST "when", BAD_CAST date);

  //Content tag
  node =xmlNewChild(root_node,0, BAD_CAST "content",0);
  node2 = xmlNewChild(node, 0, BAD_CAST "mesh",0);
  xmlNewProp(node2, BAD_CAST "name", BAD_CAST _collection->getName().c_str());

  //Splitting tag
  node=xmlNewChild(root_node,0,BAD_CAST "splitting",0);
  node2=xmlNewChild(node,0,BAD_CAST "subdomain",0);
  sprintf(buff, "%d", (int)_collection->getMesh().size());
  xmlNewProp(node2, BAD_CAST "number", BAD_CAST buff);
  node2=xmlNewChild(node,0,BAD_CAST "global_numbering",0);
  xmlNewProp(node2, BAD_CAST "present", BAD_CAST "yes");

  //Files tag
  xmlNodePtr file_node=xmlNewChild(root_node,0,BAD_CAST "files",0);

  //Mapping tag
  node = xmlNewChild(root_node,0,BAD_CAST "mapping",0);
  xmlNodePtr mesh_node = xmlNewChild(node, 0, BAD_CAST "mesh",0);
  xmlNewProp(mesh_node, BAD_CAST "name", BAD_CAST _collection->getName().c_str());

  int nbdomains= _collection->getMesh().size();

  //loop on the domains
  std::string finalMeshName="";
  if (MyGlobals::_General_Informations.size()!=0)
    {
      std::size_t found=MyGlobals::_General_Informations[0].find("finalMeshName=");
     if ((found!=std::string::npos) && (found>0))
       {
         finalMeshName=ExtractFromDescription(MyGlobals::_General_Informations[0], "finalMeshName=");
       }
    }
  if (finalMeshName.empty())
    {
      finalMeshName=_collection->getName();
    }
  for (int idomain=nbdomains-1; idomain>=0;idomain--)
    {
      std::string distfilename;
      std::ostringstream suffix;
      suffix<<filename<<idomain+1<<".med";
      distfilename=suffix.str();

      if ( !domainSelector || domainSelector->isMyDomain( idomain ) )
        {
          if ( (_collection->getMesh())[idomain]->getNumberOfCells()==0 )
            continue; //empty domain
          if (MyGlobals::_Verbose>1)
            std::cout << "proc "<< domainSelector->rank() << " : writeMedFile " << distfilename
                      << " "<< (_collection->getMesh())[idomain]->getNumberOfCells() << " cells"
                      << " " << (_collection->getFaceMesh())[idomain]->getNumberOfCells() << " faces"
                      << " " << (_collection->getMesh())[idomain]->getNumberOfNodes()<<" nodes" << std::endl;
          writeMedFile(idomain,distfilename);
        }

      if (domainSelector->rank()==0)
        {
          //updating the ascii description file
          node = xmlNewChild(file_node, 0, BAD_CAST "subfile",0);
          sprintf (buff,"%d",idomain+1);
          xmlNewProp(node, BAD_CAST "id", BAD_CAST buff);
          xmlNewChild(node,0,BAD_CAST "name",BAD_CAST distfilename.c_str());
          xmlNewChild(node,0,BAD_CAST "machine",BAD_CAST "localhost");

          node = xmlNewChild(mesh_node,0, BAD_CAST "chunk",0);
          xmlNewProp(node, BAD_CAST "subdomain", BAD_CAST buff);
          xmlNewChild(node,0,BAD_CAST "name", BAD_CAST finalMeshName.c_str());
          //xmlNewChild(node,0,BAD_CAST "name", BAD_CAST ((_collection->getMesh())[idomain]->getName()).c_str());
        }
    }
  
  //create the ascii description file
  if (domainSelector->rank()==0)
    {
      std::string myfile(filename);
      myfile.append(".xml");
      if ( !domainSelector || domainSelector->rank() == 0 )
        xmlSaveFormatFileEnc(myfile.c_str(), master_doc, "UTF-8", 1);
    }
  xmlFreeDoc(master_doc);
  xmlCleanupParser();
}
