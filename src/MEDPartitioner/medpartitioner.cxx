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

/*
  examples of launch
  export verb=1
  medpartitioner --input-file=blade.med --output-file=ttmp1_ --ndomains=2 --dump-cpu-memory --verbose=$verb
  medpartitioner --input-file=medpartitioner_blade.xml --output-file=ttmp1_ --ndomains=2 --dump-cpu-memory --verbose=$verb
  medpartitioner --input-file=ttmp1_.xml --output-file=tttmp1_ --ndomains=4 --dump-cpu-memory --verbose=$verb
*/

/*
#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_Utils.hxx"
*/

#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include <string>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace MEDPARTITIONER;

int main(int argc, char** argv)
{
#if !defined(MED_ENABLE_METIS) && !defined(MED_ENABLE_SCOTCH)
  std::cout << "Sorry, no one split method is available. Please, compile with METIS or SCOTCH." << std::endl;
  return 1;
#else

  // Defining options
  // by parsing the command line
  
  bool split_family=false;
  bool empty_groups=false;
  bool mesure_memory=false;
  bool filter_face=true;

  string input;
  string output;
  string meshname;
  string library="metis";  //default
  int ndomains;
  int help=0;
  
  //sequential : no MPI
  MyGlobals::_World_Size=1;
  MyGlobals::_Rank=0;
  MyGlobals::_Create_Boundary_Faces=0;
  MyGlobals::_Create_Joints=0;

  // Primitive parsing of command-line options
  string desc ("Available options of medpartitioner V1.0:\n"
               "\t--help                   : produces this help message\n"
               "\t--verbose                : echoes arguments\n"
               "\t--input-file=<string>    : name of the input .med file or .xml master file\n"
               "\t--output-file=<string>   : name of the resulting file (without extension)\n"
               "\t--ndomains=<number>      : number of subdomains in the output file, default is 1\n"
#if defined(MED_ENABLE_METIS) && defined(MED_ENABLE_SCOTCH)
  //user can choose!
               "\t--split-method=<string>  : name of the splitting library (metis/scotch), default is metis\n"
#endif
               "\t--create-boundary-faces  : creates boundary faces mesh in the output files\n"
               "\t--create-joints          : creates joints in the output files\n"
               "\t--dump-cpu-memory        : dumps passed CPU time and maximal increase of used memory\n"
               );

  if (argc<=1) help=1;
  string value;
  for (int i = 1; i < argc; i++)
    {
      if (strlen(argv[i]) < 3) 
        {
          cerr << "bad argument : "<< argv[i] << endl;
          return 1;
        }
    
      if (TestArg(argv[i],"--verbose",value)) 
        {
          MyGlobals::_Verbose=1;
          if (value!="") MyGlobals::_Verbose = atoi(value.c_str());
        }
      else if (TestArg(argv[i],"--help",value)) help=1;
//      else if (TestArg(argv[i],"--test",value)) test=1;
      else if (TestArg(argv[i],"--input-file",value)) input=value;
      else if (TestArg(argv[i],"--output-file",value)) output=value;
      else if (TestArg(argv[i],"--split-method",value)) library=value;
      else if (TestArg(argv[i],"--ndomains",value)) ndomains=atoi(value.c_str());
      else if (TestArg(argv[i],"--create-boundary-faces",value)) MyGlobals::_Create_Boundary_Faces=1;
      else if (TestArg(argv[i],"--create-joints",value)) MyGlobals::_Create_Joints=1;
      else if (TestArg(argv[i],"--dump-cpu-memory",value)) mesure_memory=true;
      else 
        {
          cerr << "unknown argument : "<< argv[i] << endl;
          return 1;
        }
    }

  MyGlobals::_Is0verbose=MyGlobals::_Verbose;
  
//no choice
#if defined(MED_ENABLE_METIS) && !defined(MED_ENABLE_SCOTCH)
  library = "metis";
#endif
#if !defined(MED_ENABLE_METIS) && defined(MED_ENABLE_SCOTCH)
  library = "scotch";
#endif
//user choice
#if defined(MED_ENABLE_METIS) && defined(MED_ENABLE_SCOTCH)
  if ((library!="metis") && (library!="scotch"))
    {
      cerr << "split-method only available : metis, scotch" << endl;
      return 1;
    }
#endif
 
  if (help==1)
    {
      cout<<desc<<"\n";
      return 0;
    }
  
  if (MyGlobals::_Is0verbose)
    {
      cout << "medpartitioner V1.0 :" << endl;
      cout << "  input-file = " << input << endl;
      cout << "  output-file = " << output << endl;
      cout << "  split-method = " << library << endl;
      cout << "  ndomains = " << ndomains << endl;
      cout << "  create_boundary_faces = " << MyGlobals::_Create_Boundary_Faces << endl;
      cout << "  create-joints = " << MyGlobals::_Create_Joints<< endl;
      cout << "  dump-cpu-memory = " << mesure_memory<< endl;
      cout << "  verbose = " << MyGlobals::_Verbose << endl;
    }
  //testing whether it is possible to write a file at the specified location
  if (MyGlobals::_Rank==0)
    {
      string outputtest = output + ".testioms.";
      ofstream testfile (outputtest.c_str());
      if (testfile.fail())
        { 
          cerr << "output-file directory does not exist or is in read-only access" << endl;
          return 1;
        }
      //deletes test file
      remove(outputtest.c_str());
    }
  
  // Beginning of the computation

  // Loading the mesh collection
  if (MyGlobals::_Is0verbose) cout << "Reading input files "<<endl;
  
  try
    {
      /*MEDPARTITIONER::ParaDomainSelector parallelizer(mesure_memory);
      MEDPARTITIONER::MeshCollection collection(input,parallelizer);
      MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
      aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());
      //int nbfiles=MyGlobals::_fileMedNames->size(); //nb domains
      //to have unique valid fields names/pointers/descriptions for partitionning
      collection.prepareFieldDescriptions();
      //int nbfields=collection.getFieldDescriptions().size(); //on all domains
      //cout<<ReprVectorOfString(collection.getFieldDescriptions());
    
      if (MyGlobals::_Is0verbose)
        {
          cout<<"fileNames :"<<endl
              <<ReprVectorOfString(MyGlobals::_File_Names);
          cout<<"fieldDescriptions :"<<endl
              <<ReprFieldDescriptions(collection.getFieldDescriptions()," "); //cvwat07
          cout<<"familyInfo :\n"
              <<ReprMapOfStringInt(collection.getFamilyInfo())<<endl;
          cout<<"groupInfo :\n"
              <<ReprMapOfStringVectorOfString(collection.getGroupInfo())<<endl;
        }*/
      MEDPARTITIONER::ParaDomainSelector parallelizer(mesure_memory);
      MEDPARTITIONER::MeshCollection collection(input,parallelizer);
      MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
      aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());
      //to have unique valid fields names/pointers/descriptions for partitionning
      collection.prepareFieldDescriptions();
      
      //MEDPARTITIONER::MeshCollection collection(input);

      //Creating the graph and partitioning it
      if (MyGlobals::_Is0verbose) cout << "Computing partition with " << library << endl;
  
      auto_ptr< MEDPARTITIONER::Topology > new_topo;
      if (library == "metis")
        new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::METIS));
      else
        new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH));
      parallelizer.evaluateMemory();
  
      //Creating a new mesh collection from the partitioning
      if (MyGlobals::_Is0verbose)  cout << "Creating new meshes"<< endl;
      MEDPARTITIONER::MeshCollection new_collection(collection,new_topo.get(),split_family,empty_groups);
  
      if (filter_face) new_collection.filterFaceOnCell();
    
      //to get infos on all procs
    
      //see meshName
      vector<string> finalInformations;
      vector<string> r1,r2;
      //r1=AllgathervVectorOfString(MyGlobals::_General_Informations);
      r1=MyGlobals::_General_Informations;
      //if (MyGlobals::_Is0verbose>1000) cout << "generalInformations : \n"<<ReprVectorOfString(r1);
      r2=SelectTagsInVectorOfString(r1,"ioldDomain=");
      r2=SelectTagsInVectorOfString(r2,"meshName=");
      if (r2.size()==(collection.getMesh()).size())
        {
          for (std::size_t i=0; i<r2.size(); i++)
            r2[i]=EraseTagSerialized(r2[i],"ioldDomain=");
          r2=DeleteDuplicatesInVectorOfString(r2);
          if (r2.size()==1)
            {
              string finalMesh="finalMeshName="+ExtractFromDescription(r2[0], "meshName=");
              finalInformations.push_back(SerializeFromString(finalMesh));
            }
        }
      if (finalInformations.size()==0)
        {
          if (MyGlobals::_Rank==0)
            cerr<<"Problem on final meshName : set at 'Merge'"<<endl;
          finalInformations.push_back(SerializeFromString("finalMeshName=Merge"));
        }
    
      //see field info nbComponents & componentInfo (if fields present)
      r2=SelectTagsInVectorOfString(r1,"fieldName=");
      r2=SelectTagsInVectorOfString(r2,"nbComponents=");
      //may be yes? or not?
      for (std::size_t i=0; i<r2.size(); i++)
        r2[i]=EraseTagSerialized(r2[i],"ioldFieldDouble=");
      r2=DeleteDuplicatesInVectorOfString(r2);
      for (std::size_t i=0; i<r2.size(); i++)
        finalInformations.push_back(r2[i]);
    
      MyGlobals::_General_Informations=finalInformations;
      if (MyGlobals::_Is0verbose) 
        cout << "generalInformations : \n"<<ReprVectorOfString(finalInformations);
    
      //new_collection.setSubdomainBoundaryCreates(create_boundary_faces);
      if (MyGlobals::_Is0verbose) cout << "Writing "<<ndomains<<" output files "<<output<<"xx.med"<<" and "<<output<<".xml"<<endl;
      new_collection.write(output);
  
      /*if ( mesure_memory )
        if ( parallelizer.isOnDifferentHosts() || MyGlobals::_Rank==0 )
          {
            cout << "Elapsed time = " << parallelizer.getPassedTime()
                 << ", max memory usage = " << parallelizer.evaluateMemory() << " KB"
                 << endl;
          }*/
      
      if (MyGlobals::_Is0verbose>0) cout<<"OK END"<< endl;
      return 0;
    }
  catch(const char *mess)
    {
      cerr<<mess<<endl;
      fflush(stderr);
      return 1;
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      cerr<<"INTERP_KERNEL_Exception : "<<e.what()<<endl;
      fflush(stderr);
      return 1;
    }
  catch(std::exception& e)
    {
      cerr<<"std_Exception : "<<e.what()<<endl;
      fflush(stderr);
      return 1;
    }
  catch(...)
    {
      cerr<<"an unknown type exception error was occured"<<endl;
      fflush(stderr);
      return 1;
    }
#endif
}
