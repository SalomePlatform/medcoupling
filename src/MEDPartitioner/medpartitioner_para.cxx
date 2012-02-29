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
//  Module : MED
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


#include "MEDPARTITIONER_MESHCollection.hxx"
//#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDLoader.hxx"
#include "MEDPARTITIONER_utils.hxx"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#ifdef HAVE_MPI2
#include <mpi.h>
#endif

#ifdef BOOST_PROGRAM_OPTIONS_LIB
#include <boost/program_options.hpp>
namespace po=boost::program_options;
#endif

using namespace std;
using namespace MEDPARTITIONER;

int main(int argc, char** argv)
{
#ifndef ENABLE_PARMETIS
#ifndef ENABLE_PTSCOTCH
  cout << "Sorry, no one split method is available. Please, compile with ParMETIS or PT-SCOTCH."<<endl;
  return 1;
#endif
#endif

  // Defining options
  // by parsing the command line
  //bool xml_output_master=true;
  bool split_family=false;
  bool empty_groups=false;
  bool mesure_memory=false;
  bool filter_face=true;

  string input;
  string output;
  string meshname;
  string library;
  int ndomains;
  int help=0;
  int test=0;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MyGlobals::_world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyGlobals::_rank);
  //cout<<"proc "<<MyGlobals::_rank<<" of "<<MyGlobals::_world_size<<endl; //cvw for debug
  //testVectorOfStringMPI(); //cvw
  //testRandomize();
  
  // Primitive parsing of command-line options
  string desc ("Available options of medpartitioner_para V1.0:\n"
               "\t--help                   : produces this help message\n"
               "\t--verbose                : echoes arguments\n"
               "\t--input-file=<string>    : name of the input .med file or .xml master file\n"
               "\t--output-file=<string>   : name of the resulting file (without exension)\n"
               "\t--ndomains=<number>      : number of subdomains in the output file, default is 1\n"
#ifdef ENABLE_PARMETIS
#ifdef ENABLE_PTSCOTCH
               "\t--split-method=<string>  : name of the splitting library (metis/scotch), default is metis\n"
#endif
#endif
               "\t--creates-boundary-faces : creates boundary faces mesh in the output files\n"
               //"\t--family-splitting       : preserves the family names instead of focusing on the groups\n"
               "\t--dump-cpu-memory        : dumps passed CPU time and maximal increase of used memory\n"
               //"\t--randomize=<number>     : random seed for other partitionning (only on one proc)\n"
               //"\t--atomize              : do the opposite of a good partitionner (only on one proc)\n"
               );

  if (argc<=1) help=1;
  string value;
  for (int i = 1; i < argc; i++)
  {
    if (strlen(argv[i]) < 3) 
    {
      if (MyGlobals::_rank==0) cerr << "bad argument : "<< argv[i] << endl;
      MPI_Finalize(); return 1;
    }
    
    if (testArg(argv[i],"--verbose",value)) 
    {
      MyGlobals::_verbose=1;
      if (value!="") MyGlobals::_verbose = atoi(value.c_str());
    }
    else if (testArg(argv[i],"--help",value)) help=1;
    else if (testArg(argv[i],"--test",value)) test=1;
    else if (testArg(argv[i],"--input-file",value)) input=value;
    else if (testArg(argv[i],"--output-file",value)) output=value;
    else if (testArg(argv[i],"--split-method",value)) library=value;
    //else if (testArg(argv[i],"--family-splitting",value)) split_family=true;
    else if (testArg(argv[i],"--ndomains",value)) ndomains=atoi(value.c_str());
    else if (testArg(argv[i],"--randomize",value)) MyGlobals::_randomize=atoi(value.c_str());
    else if (testArg(argv[i],"--atomize",value)) MyGlobals::_atomize=atoi(value.c_str());
    else if (testArg(argv[i],"--creates-boundary-faces",value)) MyGlobals::_creates_boundary_faces=1;
    //else if (testArg(argv[i],"--empty-groups",value)) empty_groups=true;
    else if (testArg(argv[i],"--dump-cpu-memory",value)) mesure_memory=true;
    else 
    {
      if (MyGlobals::_rank==0) cerr << "unknown argument : "<< argv[i] << endl;
      MPI_Finalize(); return 1;
    }
  }


  if (MyGlobals::_randomize!=0 && MyGlobals::_world_size!=1)
  {
    if (MyGlobals::_rank==0) cerr << "randomize only available on 1 proc (mpirun -np 1)" << endl;
    MyGlobals::_randomize=0;
  }
  
#ifdef ENABLE_PARMETIS
#ifndef ENABLE_PTSCOTCH
  library = "metis";
#endif
#else
  library = "scotch";
#endif
 
  if (help==1)
  {
    if (MyGlobals::_rank==0) cout<<desc<<"\n";
    MPI_Finalize(); return 0;
  }
  
  MyGlobals::_is0verbose=0;
  if (MyGlobals::_rank==0) MyGlobals::_is0verbose=MyGlobals::_verbose;
  //MyGlobals::_is0verbose=((MyGlobals::_rank==0) && MyGlobals::_verbose);
  if (test==1) //only for debugger
  {
    if (MyGlobals::_is0verbose>0) cout<<"tests on "<<MyGlobals::_atomize<<" "<<ndomains<<endl;
    //testPersistantMpi0To1(MyGlobals::_atomize, ndomains);
    //testPersistantMpiRing(MyGlobals::_atomize, ndomains);
    testPersistantMpiRingOnCommSplit(MyGlobals::_atomize, ndomains);
    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
    testVectorOfStringMPI();
    testMapOfStringIntMPI();
    testMapOfStringVectorOfStringMPI();
    testDataArrayMPI();
    MPI_Finalize();
    return 1;
  }
  
  if (MyGlobals::_is0verbose)
  {
    cout << "medpartitioner_para V1.0 :" << endl;
    cout << "  input-file = " << input << endl;
    cout << "  output-file = " << output << endl;
    cout << "  split-method = " << library << endl;
    //cout << "  family-splitting = " << split_family << endl;
    cout << "  ndomains = " << ndomains << endl;
    //cout << "  xml_output_master = " << xml_output_master << endl;
    cout << "  creates_boundary_faces = " << MyGlobals::_creates_boundary_faces << endl;
    //cout << "  empty_groups = " << empty_groups<< endl;
    cout << "  dump-cpu-memory = " << mesure_memory<< endl;
    cout << "  verbose = " << MyGlobals::_verbose << endl;
    //cout << "  randomize = " << MyGlobals::_randomize << endl;
    cout << "  verbose = " << MyGlobals::_verbose << endl;
  }
  //testing whether it is possible to write a file at the specified location
  if (MyGlobals::_rank==0)
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
  if (MyGlobals::_is0verbose) cout << "Reading input files "<<endl;
  
  try
  {
    MEDPARTITIONER::ParaDomainSelector parallelizer(mesure_memory);
    MEDPARTITIONER::MESHCollection collection(input,parallelizer); //cvwat01
    MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
    aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());
    //int nbfiles=MyGlobals::_fileMedNames->size(); //nb domains
    //to have unique valid fields names/pointers/descriptions for partitionning
    collection.prepareFieldDescriptions();
    //int nbfields=collection.getFieldDescriptions().size(); //on all domains
    //cout<<reprVectorOfString(collection.getFieldDescriptions());
    
    if (MyGlobals::_is0verbose)
    {
      cout<<"fileNames :"<<endl
          <<reprVectorOfString(MyGlobals::_fileNames);
      cout<<"fieldDescriptions :"<<endl
          <<reprFieldDescriptions(collection.getFieldDescriptions()," "); //cvwat07
      cout<<"familyInfo :\n"
          <<reprMapOfStringInt(collection.getFamilyInfo())<<endl;
      cout<<"groupInfo :\n"
          <<reprMapOfStringVectorOfString(collection.getGroupInfo())<<endl;
    }
    
    // Creating the graph and partitioning it   
    if (MyGlobals::_is0verbose) cout << "Computing partition "<<endl; //cvw
  
    auto_ptr< MEDPARTITIONER::Topology > new_topo;
    if (library == "metis") //cvwat06
      new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::METIS));
    else
      new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH));
    parallelizer.evaluateMemory();
  
    // Creating a new mesh collection from the partitioning
    if (MyGlobals::_is0verbose)  cout << "Creating new meshes"<<endl; //cvwat04
    MEDPARTITIONER::MESHCollection new_collection(collection,new_topo.get(),split_family,empty_groups);
    //cout<<"proc "<<MyGlobals::_rank<<" : new_collection done"<<endl;
    parallelizer.evaluateMemory();
  
    //if (!xml_output_master) new_collection.setDriverType(MEDPARTITIONER::MedAscii);
    if (filter_face) new_collection.filterFaceOnCell();
    
    //to get infos on all procs
    
    //see meshName
    vector<string> finalInformations;
    vector<string> r1,r2;
    r1=allgathervVectorOfString(MyGlobals::_generalInformations);
    //if (MyGlobals::_is0verbose>1000) cout << "generalInformations : \n"<<reprVectorOfString(r1);
    r2=selectTagsInVectorOfString(r1,"ioldDomain=");
    r2=selectTagsInVectorOfString(r2,"meshName=");
    if (r2.size()==(collection.getMesh()).size())
    {
      for (int i=0; i<r2.size(); i++) r2[i]=eraseTagSerialized(r2[i],"ioldDomain=");
      r2=deleteDuplicatesInVectorOfString(r2);
      if (r2.size()==1)
      {
        string finalMesh="finalMeshName="+extractFromDescription(r2[0], "meshName=");
        finalInformations.push_back(serializeFromString(finalMesh));
      }
    }
    if (finalInformations.size()==0)
    {
      if (MyGlobals::_rank==0)
        cerr<<"Problem on final meshName : set at 'Merge'"<<endl;
      finalInformations.push_back("finalMeshName=Merge");
    }
    
    //see field info nbComponents & componentInfo (if fields present)
    r2=selectTagsInVectorOfString(r1,"fieldName=");
    r2=selectTagsInVectorOfString(r2,"nbComponents=");
    //may be yes? or not?
    for (int i=0; i<r2.size(); i++) r2[i]=eraseTagSerialized(r2[i],"ioldFieldDouble=");
    r2=deleteDuplicatesInVectorOfString(r2);
    for (int i=0; i<r2.size(); i++) finalInformations.push_back(r2[i]);
    
    MyGlobals::_generalInformations=finalInformations;
    if (MyGlobals::_is0verbose) 
      cout << "generalInformations : \n"<<reprVectorOfString(finalInformations);
    
    //new_collection.setSubdomainBoundaryCreates(creates_boundary_faces);
    if (MyGlobals::_is0verbose) cout << "Writing "<<ndomains<<" output files "<<output<<"xx.med"<<" and "<<output<<".xml"<<endl;
    new_collection.write(output);
  
    if ( mesure_memory )
      if ( parallelizer.isOnDifferentHosts() || MyGlobals::_rank==0 )
      {
        cout << "Elapsed time = " << parallelizer.getPassedTime()
             << ", max memory usage = " << parallelizer.evaluateMemory() << " KB"
             << endl;
      }
  //   Casting the fields on the new collection
  //   if (!mesh_only)
  //     new_collection.castAllFields(*collection);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MyGlobals::_is0verbose>0) cout<<"OK END"<< endl;
    MPI_Finalize();
    return 0;
  }
  catch(const char *mess)
  {
    cerr<<"proc "<<MyGlobals::_rank<<" : "<<mess<<endl;
    fflush(stderr);
    MPI_Finalize();
    return 1;
  }
  catch(INTERP_KERNEL::Exception& e)
  {
    cerr<<"proc "<<MyGlobals::_rank<<" : INTERP_KERNEL_Exception : "<<e.what()<<endl;
    fflush(stderr);
    MPI_Finalize();
    return 1;
  }
  catch(std::exception& e)
  {
    cerr<<"proc "<<MyGlobals::_rank<<" : std_Exception : "<<e.what()<<endl;
    fflush(stderr);
    MPI_Finalize();
    return 1;
  }
  catch(...)
  {
    cerr<<"proc "<<MyGlobals::_rank<<" : an unknown type exception error was occured"<<endl;
    fflush(stderr);
    MPI_Finalize();
    return 1;
  }
}
