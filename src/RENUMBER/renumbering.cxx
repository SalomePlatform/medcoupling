// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "MEDFileData.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"

#include "MEDCouplingUMesh.hxx"

#include "RenumberingFactory.hxx"

#include <time.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;
using namespace ParaMEDMEM;
using namespace MED_RENUMBER;

int main(int argc, char** argv)
{
  double t_begin,t_read_st,t_compute_graph,t_family,t_field;
  t_begin=clock();
  if (argc <5)
    {
      cerr << "Usage : " << argv[0] 
           << " filename_in meshname method[BOOST/METIS] filename_out" << endl << endl;
      return -1;
    }
  string filename_in = argv[1];
  string meshname = argv[2];
  string type_renum = argv[3];
  string filename_out = argv[4];

  if(type_renum!="METIS" && type_renum!="BOOST")
    {
      cout << "The method has to be METIS or BOOST!" << endl;
      return -1;
    }
  // Reading file structure
  cout << "Reading : " << flush;
  MEDCouplingAutoRefCountObjectPtr<MEDFileData> fd=MEDFileData::New(filename_in.c_str());
  MEDFileMesh *m=fd->getMeshes()->getMeshWithName(meshname.c_str());
  MEDFileUMesh *mc=dynamic_cast<MEDFileUMesh *>(m);
  if(!mc)
    {
      std::ostringstream oss; oss << "In file \"" << filename_in << "\" the mesh name \"" << meshname<< "\" exists but is not unstructured !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  t_read_st=clock();
  cout << (t_read_st-t_begin)/(double) CLOCKS_PER_SEC << "s" << endl << flush;
  // Reading mesh
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> workMesh=mc->getMeshAtLevel(0);
  std::vector<int> code=workMesh->getDistributionOfTypes();
  cout << "Building the graph : " << flush;
  DataArrayInt *neighb=0,*neighbI=0;
  workMesh->computeNeighborsOfCells(neighb,neighbI);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> neighbSafe(neighb),neighbISafe(neighbI);
  const int *graph=neighbSafe->begin();
  const int *graph_index=neighbISafe->begin();
  // Compute permutation iperm->new2old perm->old2new
  vector<int> iperm,perm;
  Renumbering *renumb=RenumberingFactory(type_renum);
  renumb->renumber(graph,graph_index,workMesh->getNumberOfCells(),iperm,perm);
  delete renumb;
  t_compute_graph=clock();
  cout << " : " << (t_compute_graph-t_read_st)/(double) CLOCKS_PER_SEC << "s" << endl;
  cout.flush();
  // Connectivity
  cout << "Reordering connectivity & families and writing : " << flush;
  workMesh->renumberCells(&perm[0],false);
  mc->setMeshAtLevel(0,workMesh);
  const DataArrayInt *famField=mc->getFamilyFieldAtLevel(0);
  if(famField)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famField2=famField->renumber(&perm[0]);
      mc->setFamilyFieldArr(0,famField2);
    }
  mc->write(filename_out.c_str(),2);
  t_family=clock();
  cout << " : " << (t_family-t_compute_graph)/(double) CLOCKS_PER_SEC << "s" << endl << flush;
  // Fields
  cout << "Reordering fields and writing : " << flush;
  MEDFileFields *fs=fd->getFields();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=DataArrayInt::New();
  o2n->useArray(&perm[0],false,CPP_DEALLOC,perm.size(),1);
  fs->renumberEntitiesLyingOnMesh(meshname.c_str(),code,code,o2n);
  fs->write(filename_out.c_str(),0);
  t_field=clock();
  cout << " : " << (t_field-t_family)/(double) CLOCKS_PER_SEC << "s" << endl << flush;
  return 0;
}
