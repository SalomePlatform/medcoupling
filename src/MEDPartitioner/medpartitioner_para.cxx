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
  rm ttmp* tttmp*
  export verb=11
  mpirun -np 2 medpartitioner_para --input-file=blade.med --output-file=ttmp1_ --ndomains=2 --dump-cpu-memory --verbose=$verb
  mpirun -np 5 medpartitioner_para --input-file=blade.med --output-file=ttmp1_ --ndomains=2 --dump-cpu-memory --verbose=$verb
  mpirun -np 4 medpartitioner_para --input-file=ttmp1_.xml --output-file=tttmp1_ --ndomains=4 --verbose=$verb

  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=tmp_testMesh_20x30x50.med --output-file=ttmp2_ --verbose=$verb
  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=ttmp2_.xml --output-file=ttmp3_ --verbose=$verb

  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=tmp_testMeshWithFaces_20x30x50.med --output-file=ttmp2_ --verbose=$verb
  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=ttmp2_.xml --output-file=ttmp3_ --verbose=$verb

  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=tmp_testMesh_20x30x50_WithVecFieldOnCells.med --output-file=ttmp2_ --verbose=$verb
  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=tmp_testMesh_20x30x50_WithVecFieldOnNodes.med --output-file=ttmp2_ --verbose=$verb
  mpirun -np 2 medpartitioner_para --ndomains=2 --input-file=ttmp2_.xml --output-file=ttmp3_ --verbose=$verb

  mpirun -np 4 medpartitioner_para --ndomains=4 --input-file=tmp_testMeshHuge_20x30x50_6.xml --output-file=ttmp3_ --verbose=1
*/


#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDLoader.hxx"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;
using namespace MEDPARTITIONER;

int main(int argc, char** argv)
{
#if !defined(MED_ENABLE_PARMETIS)
  cout << "Sorry, no one split method is available. Please, compile with ParMETIS."<<endl;
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
  int test=0;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MyGlobals::_World_Size);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyGlobals::_Rank);
  
  //cout<<"proc "<<MyGlobals::_Rank<<" of "<<MyGlobals::_World_Size<<endl; //for debug
  //TestVectorOfStringMpi(); //cvw
  //TestRandomize();
  
  //Primitive parsing of command-line options
  string desc ("Available options of medpartitioner_para V1.0.3:\n"
               "\t--help                   : produces this help message\n"
               "\t--verbose                : echoes arguments\n"
               "\t--input-file=<string>    : name of the input .med file or .xml master file\n"
               "\t--output-file=<string>   : name of the resulting file (without extension)\n"
               "\t--ndomains=<number>      : number of subdomains in the output file, default is 1\n"
#if defined(MED_ENABLE_PARMETIS) 
// || defined(MED_ENABLE_PTSCOTCH)
//user can choose! (not yet)
               "\t--split-method=<string>  : name of the splitting library (metis/scotch), default is metis\n"
#endif
               "\t--create-boundary-faces : creates boundary faces mesh in the output files\n"
               "\t--dump-cpu-memory        : dumps passed CPU time and maximal increase of used memory\n"
               //"\t--randomize=<number>     : random seed for other partitionning (only on one proc)\n"
               //"\t--atomize                : do the opposite of a good partitionner (only on one proc)\n"
               );

  if (argc<=1) help=1;
  string value;
  for (int i = 1; i < argc; i++)
    {
      if (strlen(argv[i]) < 3) 
        {
          if (MyGlobals::_Rank==0) cerr << "bad argument : "<< argv[i] << endl;
          MPI_Finalize(); return 1;
        }
    
      if (TestArg(argv[i],"--verbose",value)) 
        {
          MyGlobals::_Verbose=1;
          if (value!="") MyGlobals::_Verbose = atoi(value.c_str());
        }
      else if (TestArg(argv[i],"--help",value)) help=1;
      else if (TestArg(argv[i],"--test",value)) test=1;
      else if (TestArg(argv[i],"--input-file",value)) input=value;
      else if (TestArg(argv[i],"--output-file",value)) output=value;
      else if (TestArg(argv[i],"--split-method",value)) library=value;
      else if (TestArg(argv[i],"--ndomains",value)) ndomains=atoi(value.c_str());
      else if (TestArg(argv[i],"--randomize",value)) MyGlobals::_Randomize=atoi(value.c_str());
      else if (TestArg(argv[i],"--atomize",value)) MyGlobals::_Atomize=atoi(value.c_str());
      else if (TestArg(argv[i],"--create-boundary-faces",value)) MyGlobals::_Create_Boundary_Faces=1;
      else if (TestArg(argv[i],"--dump-cpu-memory",value)) mesure_memory=true;
      else 
        {
          if (MyGlobals::_Rank==0) cerr << "unknown argument : "<< argv[i] << endl;
          MPI_Finalize(); return 1;
        }
    }


  if (MyGlobals::_Randomize!=0 && MyGlobals::_World_Size!=1)
    {
      if (MyGlobals::_Rank==0) cerr << "randomize only available on 1 proc (mpirun -np 1)" << endl;
      MyGlobals::_Randomize=0;
    }

//no choice
#if defined(MED_ENABLE_PARMETIS) && !defined(MED_ENABLE_SCOTCH)
  library = "metis";
#endif
#if !defined(MED_ENABLE_PARMETIS) && defined(MED_ENABLE_SCOTCH)
  library = "scotch";
#endif
//user choice
#if defined(MED_ENABLE_PARMETIS) && defined(MED_ENABLE_SCOTCH)
  if ((library!="metis") && (library!="scotch"))
    {
      if (MyGlobals::_Rank==0) cerr << "split-method only available : metis, scotch" << endl;
      MPI_Finalize(); return 1;
    }
#endif
 
  if (help==1)
    {
      if (MyGlobals::_Rank==0) cout<<desc<<"\n";
      MPI_Finalize(); return 0;
    }
  
  MyGlobals::_Is0verbose=0;
  if (MyGlobals::_Rank==0) MyGlobals::_Is0verbose=MyGlobals::_Verbose;
  //MyGlobals::_Is0verbose=((MyGlobals::_Rank==0) && MyGlobals::_Verbose);
  if (test==1) //only for debugger
    {
      if (MyGlobals::_Is0verbose>0) cout<<"tests on "<<MyGlobals::_Atomize<<" "<<ndomains<<endl;
      //TestPersistantMpi0To1(MyGlobals::_Atomize, ndomains);
      //TestPersistantMpiRing(MyGlobals::_Atomize, ndomains);
      TestPersistantMpiRingOnCommSplit(MyGlobals::_Atomize, ndomains);
      //MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      return 0;
      TestVectorOfStringMpi();
      TestMapOfStringIntMpi();
      TestMapOfStringVectorOfStringMpi();
      TestDataArrayMpi();
      MPI_Finalize();
      return 1;
    }
  
  if (MyGlobals::_Is0verbose)
    {
      cout << "medpartitioner_para V1.0 :" << endl;
      cout << "  input-file = " << input << endl;
      cout << "  output-file = " << output << endl;
      cout << "  split-method = " << library << endl;
      cout << "  ndomains = " << ndomains << endl;
      cout << "  create_boundary_faces = " << MyGlobals::_Create_Boundary_Faces << endl;
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
          MPI_Finalize(); return 1;
        }
      //deletes test file
      remove(outputtest.c_str());
    }
  
  // Beginning of the computation

  // Loading the mesh collection
  if (MyGlobals::_Is0verbose) cout << "Reading input files "<<endl;
  
  try
    {
      MEDPARTITIONER::ParaDomainSelector parallelizer(mesure_memory);
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
              <<ReprFieldDescriptions(collection.getFieldDescriptions()," ");
          cout<<"familyInfo :\n"
              <<ReprMapOfStringInt(collection.getFamilyInfo())<<endl;
          cout<<"groupInfo :\n"
              <<ReprMapOfStringVectorOfString(collection.getGroupInfo())<<endl;
        }
    
      //Creating the graph and partitioning it
      if (MyGlobals::_Is0verbose) cout << "Computing partition "<< endl;
  
      auto_ptr< MEDPARTITIONER::Topology > new_topo;
      if (library == "metis")
        new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::METIS));
      else
        new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH));
      parallelizer.evaluateMemory();
  
      //Creating a new mesh collection from the partitioning
      if (MyGlobals::_Is0verbose)  cout << "Creating new meshes"<< endl;
      MEDPARTITIONER::MeshCollection new_collection(collection,new_topo.get(),split_family,empty_groups);
      //cout<<"proc "<<MyGlobals::_Rank<<" : new_collection done"<<endl;
      parallelizer.evaluateMemory();
  
      //if (!xml_output_master) new_collection.setDriverType(MEDPARTITIONER::MedAscii);
      if (filter_face) new_collection.filterFaceOnCell();
    
      //to get infos on all procs
    
      //see meshName
      vector<string> finalInformations;
      vector<string> r1,r2;
      r1=AllgathervVectorOfString(MyGlobals::_General_Informations);
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
  
      if ( mesure_memory )
        if ( parallelizer.isOnDifferentHosts() || MyGlobals::_Rank==0 )
          {
            cout << "Elapsed time = " << parallelizer.getPassedTime()
                 << ", max memory usage = " << parallelizer.evaluateMemory() << " KB"
                 << endl;
          }
      if(MyGlobals::_World_Size>1)
        MPI_Barrier(MPI_COMM_WORLD);
      if (MyGlobals::_Is0verbose>0) cout<<"OK END"<< endl;
      MPI_Finalize();
      return 0;
    }
  catch(const char *mess)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : "<<mess<<endl;
      fflush(stderr);
      MPI_Finalize();
      return 1;
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : INTERP_KERNEL_Exception : "<<e.what()<<endl;
      fflush(stderr);
      MPI_Finalize();
      return 1;
    }
  catch(std::exception& e)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : std_Exception : "<<e.what()<<endl;
      fflush(stderr);
      MPI_Finalize();
      return 1;
    }
  catch(...)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : an unknown type exception error was occured"<<endl;
      fflush(stderr);
      MPI_Finalize();
      return 1;
    }
#endif
}

