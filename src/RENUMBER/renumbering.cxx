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
using namespace MEDCoupling;
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
  MCAuto<MEDFileData> fd(MEDFileData::New(filename_in));
  MEDFileMesh *m=fd->getMeshes()->getMeshWithName(meshname);
  MEDFileUMesh *mc=dynamic_cast<MEDFileUMesh *>(m);
  if(!mc)
    {
      std::ostringstream oss; oss << "In file \"" << filename_in << "\" the mesh name \"" << meshname<< "\" exists but is not unstructured !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  t_read_st=clock();
  cout << (t_read_st-t_begin)/(double) CLOCKS_PER_SEC << "s" << endl << flush;
  // Reading mesh
  MCAuto<MEDCouplingUMesh> workMesh=mc->getMeshAtLevel(0);
  std::vector<int> code=workMesh->getDistributionOfTypes();
  cout << "Building the graph : " << flush;
  DataArrayInt *neighb=0,*neighbI=0;
  workMesh->computeNeighborsOfCells(neighb,neighbI);
  MCAuto<DataArrayInt> neighbSafe(neighb),neighbISafe(neighbI),ipermSafe,permSafe;
  const int *graph=neighbSafe->begin();
  const int *graph_index=neighbISafe->begin();
  // Compute permutation iperm->new2old perm->old2new
  DataArrayInt *iperm(0),*perm(0);
  Renumbering *renumb=RenumberingFactory(type_renum);
  renumb->renumber(graph,graph_index,workMesh->getNumberOfCells(),iperm,perm);
  ipermSafe=iperm; permSafe=perm;
  delete renumb;
  ipermSafe=0;//erase new2old, we are using only old 2 new
  t_compute_graph=clock();
  cout << " : " << (t_compute_graph-t_read_st)/(double) CLOCKS_PER_SEC << "s" << endl;
  cout.flush();
  // Connectivity
  cout << "Reordering connectivity & families and writing : " << flush;
  workMesh->renumberCells(perm->begin(),false);
  mc->setMeshAtLevel(0,workMesh);
  const DataArrayInt *famField=mc->getFamilyFieldAtLevel(0);
  if(famField)
    {
      MCAuto<DataArrayInt> famField2=famField->renumber(perm->begin());
      mc->setFamilyFieldArr(0,famField2);
    }
  mc->write(filename_out,2);
  t_family=clock();
  cout << " : " << (t_family-t_compute_graph)/(double) CLOCKS_PER_SEC << "s" << endl << flush;
  // Fields
  cout << "Reordering fields and writing : " << flush;
  MEDFileFields *fs=fd->getFields();
  if(fs)
    {
      for(int i=0;i<fs->getNumberOfFields();i++)
        {
          MEDFileFieldMultiTS *fmts=dynamic_cast<MEDFileFieldMultiTS *>(fs->getFieldAtPos(i));
          if(!fmts) continue;
          if(fmts->getMeshName()==meshname)
            {
              for(int j=0;j<fmts->getNumberOfTS();j++)
                {
                  MEDFileField1TS *f1ts=dynamic_cast<MEDFileField1TS*>(fmts->getTimeStepAtPos(j));
                  if(!f1ts) continue;
                  DataArrayDouble *arr=f1ts->getUndergroundDataArray();
                  arr->renumberInPlace(perm->begin());
                }
            }
        }
      fs->write(filename_out,0);
      //fs->renumberEntitiesLyingOnMesh(meshname,code,code,o2n); bugged
    }
  t_field=clock();
  cout << " : " << (t_field-t_family)/(double) CLOCKS_PER_SEC << "s" << endl << flush;
  return 0;
}
