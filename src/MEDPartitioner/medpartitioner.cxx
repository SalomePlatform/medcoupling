//  Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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
//  MED medsplitter : tool to split n MED files into p separate 
//                    MED files with a partitioning specified 
//                    by an external tool
//  File   : medsplitter.cxx
//  Author : Vincent BERGEAUD (CEA-DEN/DANS/DM2S/SFME/LGLS)
//  Module : MED
//

#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"

#ifdef BOOST_PROGRAM_OPTIONS_LIB
#include <boost/program_options.hpp>
#endif

#include <string>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <iostream>

#ifdef BOOST_PROGRAM_OPTIONS_LIB
namespace po=boost::program_options;
#endif

int main(int argc, char** argv)
{
#ifndef ENABLE_METIS
#ifndef ENABLE_SCOTCH
  std::cout << "Sorry, no one split method is available. Please, compile with METIS or SCOTCH." << std::endl;
  return 1;
#endif
#endif

  // Defining options
  // by parsing the command line
  bool mesh_only = false;
  bool is_sequential = true;
  bool xml_output_master=true;
  bool creates_boundary_faces=false;  
  bool split_families=false;
  bool empty_groups=false;

  std::string input;
  std::string output;
  std::string meshname;
  std::string library;
  int ndomains;

#ifdef BOOST_PROGRAM_OPTIONS_LIB
  // Use boost::program_options for command-line options parsing
  po::options_description desc("Available options of medpartitioner V1.0");
  desc.add_options()
    ("help","produces this help message")
    ("mesh-only","prevents the splitter from creating the fields contained in the original file(s)")
    ("distributed","specifies that the input file is distributed")
    ("input-file",po::value<std::string>(),"name of the input MED file")
    ("output-file",po::value<std::string>(),"name of the resulting file")
    ("meshname",po::value<std::string>(),"name of the input mesh")
#ifdef ENABLE_METIS
#ifdef ENABLE_SCOTCH
    ("split-method",po::value<std::string>(&library)->default_value("metis"),"name of the splitting library (metis,scotch)")
#endif
#endif
    ("ndomains",po::value<int>(&ndomains)->default_value(1),"number of subdomains in the output file")
    ("plain-master","creates a plain masterfile instead of an XML file")
    ("creates-boundary-faces","creates the necessary faces so that faces joints are created in the output files")
    ("family-splitting","preserves the family names instead of focusing on the groups")
    ("empty-groups","creates empty groups in zones that do not contain a group from the original domain");

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);

  if (vm.count("help"))
    {
      std::cout << desc << "\n";
      return 1;
    }

  if (!vm.count("ndomains"))
    {
      std::cout << "ndomains must be specified !"<< std::endl;
      return 1;
    }

  ndomains = vm["ndomains"].as<int>();
  if (!vm.count("input-file") || !vm.count("output-file"))
    {
      std::cout << "input-file and output-file names must be specified" << std::endl;
      return 1;
    }

  if (!vm.count("distributed") && !vm.count("meshname") )
    {
      std::cout << "for a serial MED file, mesh name must be selected with --meshname=..." << std::endl;
      return 1;
    }

  input = vm["input-file"].as<std::string>();
  output = vm["output-file"].as<std::string>();

  if (vm.count("mesh-only"))
    mesh_only=true;

  if (vm.count("distributed"))
    is_sequential=false;

  if (is_sequential)
    meshname = vm["meshname"].as<std::string>();

  if (vm.count("plain-master"))
    xml_output_master=false;

  if (vm.count("creates-boundary-faces"))
    creates_boundary_faces=true;

  if (vm.count("split-families"))
    split_families=true;

  if (vm.count("empty-groups"))
    empty_groups=true;

#else // BOOST_PROGRAM_OPTIONS_LIB

  // Primitive parsing of command-line options

  std::string desc ("Available options of medpartitioner V1.0:\n"
                    "\t--help                 : produces this help message\n"
                    "\t--mesh-only            : do not create the fields contained in the original file(s)\n"
                    "\t--distributed          : specifies that the input file is distributed\n"
                    "\t--input-file=<string>  : name of the input MED file\n"
                    "\t--output-file=<string> : name of the resulting file\n"
                    "\t--meshname=<string>    : name of the input mesh (not used with --distributed option)\n"
                    "\t--ndomains=<number>    : number of subdomains in the output file, default is 1\n"
#ifdef ENABLE_METIS
#ifdef ENABLE_SCOTCH
               "\t--split-method=<string>: name of the splitting library (metis/scotch), default is metis\n"
#endif
#endif
               "\t--plain-master         : creates a plain masterfile instead of an XML file\n"
               "\t--creates-boundary-faces: creates the necessary faces so that faces joints are created in the output files\n"
               "\t--family-splitting     : preserves the family names instead of focusing on the groups\n"
               "\t--empty-groups         : creates empty groups in zones that do not contain a group from the original domain"
               );

  if (argc < 4)
    {
      std::cout << desc.c_str() << std::endl;
      return 1;
    }

  for (int i = 1; i < argc; i++)
    {
      if (strlen(argv[i]) < 3)
        {
          std::cout << desc.c_str() << std::endl;
          return 1;
        }
      
      if (strncmp(argv[i],"--m",3) == 0)
        {
          if (strcmp(argv[i],"--mesh-only") == 0)
            {
              mesh_only = true;
              std::cout << "\tmesh-only = " << mesh_only << std::endl; // tmp
            }
          else if (strlen(argv[i]) > 11)
            { // "--meshname="
              meshname = (argv[i] + 11);
              std::cout << "\tmeshname = " << meshname << std::endl; // tmp
            }
        }
      else if (strncmp(argv[i],"--d",3) == 0)
        {
          is_sequential = false;
          std::cout << "\tis_sequential = " << is_sequential << std::endl; // tmp
        }
      else if (strncmp(argv[i],"--i",3) == 0)
        {
          if (strlen(argv[i]) > 13)
            { // "--input-file="
              input = (argv[i] + 13);
              std::cout << "\tinput-file = " << input << std::endl; // tmp
            }
        }
      else if (strncmp(argv[i],"--o",3) == 0)
        {
          if (strlen(argv[i]) > 14)
            { // "--output-file="
              output = (argv[i] + 14);
              std::cout << "\toutput-file = " << output << std::endl; // tmp
            }
        }
      else if (strncmp(argv[i],"--s",3) == 0)
        {
          if (strlen(argv[i]) > 15)
            { // "--split-method="
              library = (argv[i] + 15);
              std::cout << "\tsplit-method = " << library << std::endl; // tmp
            }
        }
      else if (strncmp(argv[i],"--f",3) == 0)
        { //"--family-splitting"
          split_families=true;
          std::cout << "\tfamily-splitting true" << std::endl; // tmp
        }
      else if (strncmp(argv[i],"--n",3) == 0)
        {
          if (strlen(argv[i]) > 11)
            { // "--ndomains="
              ndomains = atoi(argv[i] + 11);
              std::cout << "\tndomains = " << ndomains << std::endl; // tmp
            }
        }
      else if (strncmp(argv[i],"--p",3) == 0)
        { // "--plain-master"
          xml_output_master = false;
          std::cout << "\txml_output_master = " << xml_output_master << std::endl; // tmp
        }
      else if (strncmp(argv[i],"--c",3) == 0)
        { // "--creates-boundary-faces"
          creates_boundary_faces = true;
          std::cout << "\tcreates_boundary_faces = " << creates_boundary_faces << std::endl; // tmp
        }
    else if (strncmp(argv[i],"--e",3) == 0)
      { // "--empty-groups"
        empty_groups = true;
        std::cout << "\tempty_groups = true" << std::endl; // tmp
      }
    else
      {
        std::cout << desc.c_str() << std::endl;
        return 1;
      }
    }
  
  if (is_sequential && meshname.empty())
    {
      std::cout << "Mesh name must be given for sequential(not distributed) input file." << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }
  
#endif // BOOST_PROGRAM_OPTIONS_LIB


  //testing whether it is possible to write a file at the specified location
  std::string outputtest = output + ".testioms.";
  std::ofstream testfile (outputtest.c_str());
  if (testfile.fail())
    { 
      std::cout << "MEDPARTITIONER : output-file directory does not exist or is in read-only access" << std::endl;
      return 1;
    };
  //deletes test file
  remove(outputtest.c_str());

  // Beginning of the computation

  // Loading the mesh collection
  MEDPARTITIONER::MeshCollection* collection;
  std::cout << "MEDPARTITIONER : reading input files "<< std::endl;
  if (is_sequential)
    collection = new MEDPARTITIONER::MeshCollection(input,meshname);
  else
    collection = new MEDPARTITIONER::MeshCollection(input);

  std::cout << "MEDPARTITIONER : computing partition "<< std::endl;

  // Creating the graph and partitioning it   
#ifdef ENABLE_METIS
#ifndef ENABLE_SCOTCH
  library = "metis";
#endif
#else
  library = "scotch";
#endif
  std::cout << "\tsplit-method = " << library << std::endl; // tmp

  MEDPARTITIONER::Topology* new_topo;
  if (library == "metis")
    new_topo = collection->createPartition(ndomains,MEDPARTITIONER::Graph::METIS);
  else
    new_topo = collection->createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH);

  std::cout << "MEDPARTITIONER : creating new meshes"<< std::endl;

  // Creating a new mesh collection from the partitioning
  MEDPARTITIONER::MeshCollection new_collection(*collection, new_topo, split_families, empty_groups);
  if (mesh_only)
    {
      delete collection;
      collection=0;
    }

  if (!xml_output_master)
    new_collection.setDriverType(MEDPARTITIONER::MedAscii);

  //  new_collection.setSubdomainBoundaryCreates(creates_boundary_faces);

  std::cout << "MEDPARTITIONER : writing output files "<< std::endl;
  new_collection.write(output);

  // Casting the fields on the new collection
  //   if (!mesh_only)
  //     new_collection.castAllFields(*collection);


  // Cleaning memory
  delete collection;
  delete new_topo;

  return 0;
}
