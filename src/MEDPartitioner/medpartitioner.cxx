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
//  MED medsplitter : tool to split n MED files into p separate 
//                    MED files with a partitioning specified 
//                    by an external tool
//  File   : medsplitter.cxx
//  Author : Vincent BERGEAUD (CEA-DEN/DANS/DM2S/SFME/LGLS)
//  Module : MED
//
#ifdef BOOST_PROGRAM_OPTIONS_LIB
#include <boost/program_options.hpp>
namespace po=boost::program_options;
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_MESHCollection.hxx"
#include "MEDPARTITIONER_Topology.hxx"

using namespace std;

int main(int argc, char** argv)
{
#ifndef ENABLE_METIS
#ifndef ENABLE_SCOTCH
  cout << "Sorry, no one split method is available. Please, compile with METIS or SCOTCH."<<endl;
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

  string input;
  string output;
  string meshname;
  string library;
  int ndomains;

#ifdef BOOST_PROGRAM_OPTIONS_LIB
  // Use boost::program_options for command-line options parsing
  po::options_description desc("Available options of medpartitioner V1.0");
  desc.add_options()
    ("help","produces this help message")
    ("mesh-only","prevents the splitter from creating the fields contained in the original file(s)")
    ("distributed","specifies that the input file is distributed")
    ("input-file",po::value<string>(),"name of the input MED file")
    ("output-file",po::value<string>(),"name of the resulting file")
    ("meshname",po::value<string>(),"name of the input mesh")
#ifdef ENABLE_METIS
#ifdef ENABLE_SCOTCH
    ("split-method",po::value<string>(&library)->default_value("metis"),"name of the splitting library (metis,scotch)")
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
    cout<<desc<<"\n";
    return 1;
  }

  if (!vm.count("ndomains"))
  {
    cout << "ndomains must be specified !"<<endl;
    return 1;
  }

  ndomains = vm["ndomains"].as<int>();
  if (!vm.count("input-file") || !vm.count("output-file"))
  {
    cout << "input-file and output-file names must be specified"<<endl;
    return 1;
  }

  if (!vm.count("distributed") && !vm.count("meshname") )
  {
    cout << "for a serial MED file, mesh name must be selected with --meshname=..."<<endl;
    return 1;
  }

  input = vm["input-file"].as<string>();
  output = vm["output-file"].as<string>();

  if (vm.count("mesh-only"))
    mesh_only=true;

  if (vm.count("distributed"))
    is_sequential=false;

  if (is_sequential)
    meshname = vm["meshname"].as<string>();

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

  string desc ("Available options of medpartitioner V1.0:\n"
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

  if (argc < 4) {
    cout << desc.c_str() << endl;
    return 1;
  }

  for (int i = 1; i < argc; i++) {
    if (strlen(argv[i]) < 3) {
      cout << desc.c_str() << endl;
      return 1;
    }

    if (strncmp(argv[i],"--m",3) == 0) {
      if (strcmp(argv[i],"--mesh-only") == 0) {
        mesh_only = true;
        cout << "\tmesh-only = " << mesh_only << endl; // tmp
      }
      else if (strlen(argv[i]) > 11) { // "--meshname="
        meshname = (argv[i] + 11);
        cout << "\tmeshname = " << meshname << endl; // tmp
      }
    }
    else if (strncmp(argv[i],"--d",3) == 0) {
      is_sequential = false;
      cout << "\tis_sequential = " << is_sequential << endl; // tmp
    }
    else if (strncmp(argv[i],"--i",3) == 0) {
      if (strlen(argv[i]) > 13) { // "--input-file="
        input = (argv[i] + 13);
        cout << "\tinput-file = " << input << endl; // tmp
      }
    }
    else if (strncmp(argv[i],"--o",3) == 0) {
      if (strlen(argv[i]) > 14) { // "--output-file="
        output = (argv[i] + 14);
        cout << "\toutput-file = " << output << endl; // tmp
      }
    }
    else if (strncmp(argv[i],"--s",3) == 0) {
      if (strlen(argv[i]) > 15) { // "--split-method="
        library = (argv[i] + 15);
        cout << "\tsplit-method = " << library << endl; // tmp
      }
    }
    else if (strncmp(argv[i],"--f",3) == 0) { //"--family-splitting"
      split_families=true;
      cout << "\tfamily-splitting true" << endl; // tmp
    }
    else if (strncmp(argv[i],"--n",3) == 0) {
      if (strlen(argv[i]) > 11) { // "--ndomains="
        ndomains = atoi(argv[i] + 11);
        cout << "\tndomains = " << ndomains << endl; // tmp
      }
    }
    else if (strncmp(argv[i],"--p",3) == 0) { // "--plain-master"
      xml_output_master = false;
      cout << "\txml_output_master = " << xml_output_master << endl; // tmp
    }
    else if (strncmp(argv[i],"--c",3) == 0) { // "--creates-boundary-faces"
      creates_boundary_faces = true;
      cout << "\tcreates_boundary_faces = " << creates_boundary_faces << endl; // tmp
    }
    else if (strncmp(argv[i],"--e",3) == 0) { // "--empty-groups"
      empty_groups = true;
      cout << "\tempty_groups = true" << endl; // tmp
    }
    else {
      cout << desc.c_str() << endl;
      return 1;
    }
  }

  if (is_sequential && meshname.empty()) {
    cout << "Mesh name must be given for sequential(not distributed) input file." << endl;
    cout << desc << endl;
    return 1;
  }

#endif // BOOST_PROGRAM_OPTIONS_LIB


  //testing whether it is possible to write a file at the specified location
  string outputtest = output + ".testioms.";
  ofstream testfile (outputtest.c_str());
  if (testfile.fail())
  { 
    cout << "MEDPARTITIONER : output-file directory does not exist or is in read-only access" << endl;
    return 1;
  };
  //deletes test file
  remove(outputtest.c_str());

  // Beginning of the computation

  // Loading the mesh collection
  MEDPARTITIONER::MESHCollection* collection;
  cout << "MEDPARTITIONER : reading input files "<<endl;
  if (is_sequential)
    collection = new MEDPARTITIONER::MESHCollection(input,meshname);
  else
    collection = new MEDPARTITIONER::MESHCollection(input);

  cout << "MEDPARTITIONER : computing partition "<<endl;

  // Creating the graph and partitioning it   
#ifdef ENABLE_METIS
#ifndef ENABLE_SCOTCH
  library = "metis";
#endif
#else
  library = "scotch";
#endif
  cout << "\tsplit-method = " << library << endl; // tmp

  MEDPARTITIONER::Topology* new_topo;
  if (library == "metis")
    new_topo = collection->createPartition(ndomains,MEDPARTITIONER::Graph::METIS);
  else
    new_topo = collection->createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH);

  cout << "MEDPARTITIONER : creating new meshes"<<endl;

  // Creating a new mesh collection from the partitioning
  MEDPARTITIONER::MESHCollection new_collection(*collection, new_topo, split_families, empty_groups);
  if (mesh_only)
  {
    delete collection;
    collection=0;
  }

  if (!xml_output_master)
    new_collection.setDriverType(MEDPARTITIONER::MedAscii);

  //  new_collection.setSubdomainBoundaryCreates(creates_boundary_faces);

  cout << "MEDPARTITIONER : writing output files "<<endl;
  new_collection.write(output);

  // Casting the fields on the new collection
//   if (!mesh_only)
//     new_collection.castAllFields(*collection);


  // Cleaning memory
  delete collection;
  delete new_topo;

  return 0;
}
