// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "MEDMEM_Family.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Meshing.hxx"
#include "MEDMEM_MedMeshDriver.hxx"
#include "MEDMEM_Connectivity.hxx"
#include "MEDMEM_Field.hxx"
#include "MEDMEM_DriversDef.hxx"
#include "MEDMEM_MedFileBrowser.hxx"
#include "MEDMEM_MedMeshDriver.hxx"

#include "RenumberingFactory.hxx"

#include <time.h>
using namespace MEDMEM;
using namespace std;
using namespace MED_EN;
using namespace MED_RENUMBER;

void computeNeighbour(const MESH* mesh,const medGeometryElement& Type, vector<list<int> >& neighbour, int& ntot,int& nb_cell)
{
  CONNECTIVITY* conn = (CONNECTIVITY*)mesh->getConnectivityptr();
  conn->calculateFullDescendingConnectivity(MED_CELL);
  const int* rev_conn=mesh->getReverseConnectivity(MED_EN::MED_DESCENDING, MED_EN::MED_CELL);
  const int* rev_conn_index=mesh->getReverseConnectivityIndex(MED_EN::MED_DESCENDING, MED_EN::MED_CELL);
  int nb_face= mesh->getNumberOfElements(MED_FACE,MED_ALL_ELEMENTS);
  int nb_edge = mesh->getNumberOfElements(MED_EDGE,MED_ALL_ELEMENTS);
  nb_cell= mesh->getNumberOfElements(MED_CELL,Type);

  int nb_constituent;
  if(mesh->getMeshDimension()==2)
    nb_constituent = nb_edge;
  else if (mesh->getMeshDimension()==3)
    nb_constituent = nb_face;
  else
    throw MEDEXCEPTION("Wrong dimension");

  neighbour.resize(nb_cell,(list<int>)0);
  ntot=0;
  for(int i=0;i<nb_constituent;++i)
    {
      for(int j=rev_conn_index[i]-1;j<rev_conn_index[i+1]-1;++j)
        {
          for(int k=j+1;k<rev_conn_index[i+1]-1;++k)
            {
              if(rev_conn[j]!=0 && rev_conn[k]!=0)
                {
                  ntot+=2;
                  neighbour[rev_conn[j]-1].push_back(rev_conn[k]);
                  neighbour[rev_conn[k]-1].push_back(rev_conn[j]);
                }
            }
        }
    }
}

void changeConnectivity(MESH& mesh, const medGeometryElement& Type, const int& nb_cell, const vector<int>& iperm)
{
  /*if(Type==MED_POLYHEDRA)
    {
      int *conn_face_index_init=(int*)mesh.getPolyhedronFacesIndex();
      int *conn_index_init=(int*)mesh.getPolyhedronIndex(MED_FULL_INTERLACE);
      int *conn_init=(int*)mesh.getPolyhedronConnectivity(MED_FULL_INTERLACE);

      int *conn_index_renum=new int[nb_cell+1];
      int *conn_face_index_renum=new int[conn_index_init[nb_cell]];
      int *conn_renum=new int[conn_face_index_init[conn_index_init[nb_cell]-1]-1];

      int i_cell,i_face,i_conn;
      int iter_face=0;
      int iter_conn=0;
      int i2;
      conn_index_renum[0]=1;
      conn_face_index_renum[0]=1;
      for(i_cell=0;i_cell<nb_cell;++i_cell)
        {
          i2=iperm[i_cell]-1;
          for(i_face=conn_index_init[i2]-1;i_face<conn_index_init[i2+1]-1;++i_face)
            {
              for(i_conn=conn_face_index_init[i_face]-1;i_conn<conn_face_index_init[i_face+1]-1;++i_conn)
                {
                  conn_renum[iter_conn]=conn_init[i_conn];
                  ++iter_conn;
                }
              conn_face_index_renum[iter_face+1]=iter_conn+1;
              ++iter_face;
            }
          conn_index_renum[i_cell+1]=iter_face+1;
        }
      memcpy(conn_face_index_init,conn_face_index_renum,sizeof(int)*conn_index_init[nb_cell]);
      memcpy(conn_index_init,conn_index_renum,sizeof(int)*(nb_cell+1));
      memcpy(conn_init,conn_renum, sizeof(int)*(conn_face_index_init[conn_index_init[nb_cell]-1]-1));

      delete[] conn_index_renum;
      delete[] conn_face_index_renum;
      delete[] conn_renum;
    }
  else if (Type==MED_POLYGON)
    {
      int *conn_init=(int*)mesh.getPolygonsConnectivity(MED_FULL_INTERLACE,MED_CELL);
      int *conn_index_init=(int*)mesh.getPolygonsConnectivityIndex(MED_FULL_INTERLACE,MED_CELL);
      int *conn_index_renum=new int[nb_cell+1];
      int *conn_renum=new int[conn_index_init[nb_cell]-1];

      int iter=0;
      int i2;
      conn_index_renum[0]=1;
      for(int i=0;i<nb_cell;++i)
        {
          i2=iperm[i]-1;
          for(int k=conn_index_init[i2];k<conn_index_init[i2+1];++k)
            {
              conn_renum[iter]=conn_init[k-1];
              ++iter;
            }
          conn_index_renum[i+1]=iter+1;
        }
      memcpy(conn_index_init,conn_index_renum,sizeof(int)*(nb_cell+1));
      memcpy(conn_init,conn_renum, sizeof(int)*(conn_index_init[nb_cell]-1));

      delete[] conn_renum;
      delete[] conn_index_renum;
    }
    else*/
    {
      const int *conn_init=mesh.getConnectivity(MED_NODAL,MED_CELL,Type);
      const int *conn_index_init=mesh.getConnectivityIndex(MED_NODAL,MED_CELL);
      int *conn_renum=new int[conn_index_init[nb_cell]-1];
      int *conn_index_renum=new int[nb_cell+1];

      int iter=0;
      int i2;
      conn_index_renum[0]=1;
      for(int i=0;i<nb_cell;++i)
        {
          i2=iperm[i]-1;
          for(int k=conn_index_init[i2];k<conn_index_init[i2+1];++k)
            {
              conn_renum[iter]=conn_init[k-1];
              ++iter;
            }
          conn_index_renum[i+1]=iter+1;
        }

      CONNECTIVITY* myConnectivity=(CONNECTIVITY*)mesh.getConnectivityptr();
      myConnectivity->setNodal(conn_renum,MED_CELL,Type,conn_index_renum);
      delete[] conn_renum;
      delete[] conn_index_renum;
    }
}

void changeFamily(MESH* mesh, const medGeometryElement& Type, const vector<int>& perm)
{
  int nb_families=mesh->getNumberOfFamilies(MED_CELL);
  for (int i=0;i<nb_families;++i)
    {
      const FAMILY* family=mesh->getFamily(MED_CELL,i+1);
      if (!family->isOnAllElements())
        {
          int nb_elem=family->getNumberOfElements(Type);
          int *number=(int *)family->getNumber(Type);
          for(int j=0;j<nb_elem;++j)
            number[j]=perm[number[j]-1];
        }
    }
}

int main (int argc, char** argv)
{
  double t_begin,t_read_st,t_read_mesh,t_compute_graph,t_connectiv,t_family,t_field;
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
      exit(-1);
    }

  string s="rm "+filename_out;
  system(s.c_str());

  // Reading file structure
  const MEDFILEBROWSER med_struct(filename_in);
  int nb_mesh, nb_fields;
  vector<string> mesh_names,f_names;
  nb_mesh=med_struct.getNumberOfMeshes();
  nb_fields=med_struct.getNumberOfFields();
  mesh_names=med_struct.getMeshNames();
  f_names=med_struct.getFieldNames();
  if(nb_mesh!=1)
    {
      cout << "There are many meshes in the file" << endl;
      return -1;
    }
  if(mesh_names[0].c_str()!=meshname)
    {
      cout << "Mesh name does not match" << endl;
      return -1;
    }
  vector<string> field_names;
  vector<int> iternumber;
  vector<int> ordernumber;
  vector<int> types;
  int nb_fields_tot=0;
  for (int ifield = 0; ifield < nb_fields; ifield++)
    {
      vector<DT_IT_> dtit=med_struct.getFieldIteration(f_names[ifield]);
      for (vector<DT_IT_>::const_iterator iter =dtit.begin(); iter!=dtit.end(); iter++)
        {
          field_names.push_back(f_names[ifield]);
          iternumber.push_back(iter->dt);
          ordernumber.push_back(iter->it);
          ++nb_fields_tot;
          if(med_struct.getFieldType(f_names[ifield])==MED_EN::MED_REEL64)
            types.push_back(1);
          else
            types.push_back(0);

        }
    }
  t_read_st=clock();

  // Reading mesh
  MESH myMesh;
  myMesh.setName(meshname);
  MED_MESH_RDONLY_DRIVER *drv22=new MED_MESH_RDONLY_DRIVER(filename_in,&myMesh);
  drv22->desactivateFacesComputation();
  int newDrv=myMesh.addDriver(*drv22);
  delete drv22;
  myMesh.read(newDrv);
  int nb_type=myMesh.getNumberOfTypes(MED_CELL);
  if (nb_type!=1)
    {
      cout << "Mesh must have only one type of cell" << endl;
      return -1;
    }
  const medGeometryElement *Types = myMesh.getTypes(MED_CELL);
  medGeometryElement Type=Types[0];

  t_read_mesh=clock();
  MESH* workMesh=new MESH(myMesh);
  cout << "Building the graph    ";
  cout.flush();
  int ntot,nb_cell;
  vector<list<int> > neighbour;
  computeNeighbour(workMesh,Type,neighbour,ntot,nb_cell);
  int* graph=new int[ntot];
  int* graph_index=new int[nb_cell+1];
  graph_index[0]=1;
  int count=0;
  for(int i=0;i<nb_cell;++i)
    {
      for (list<int>::const_iterator it=neighbour[i].begin();it!=neighbour[i].end();++it)
        {
          graph[count]=*it;
          ++count;
        }
      graph_index[i+1]=count+1;
    }


  // Compute permutation
  vector<int> iperm,perm;
  Renumbering* renumb= RenumberingFactory(type_renum);
  renumb->renumber(graph,graph_index,nb_cell,iperm,perm);
  delete renumb;
  delete workMesh;
  t_compute_graph=clock();
  cout << " : " << (t_compute_graph-t_read_mesh)/(double) CLOCKS_PER_SEC << "s" << endl;
  cout.flush();

  // Connectivity
  cout << "Computing connectivity";
  cout.flush();
  MESH meshRenum(myMesh);
  changeConnectivity(meshRenum,Type,nb_cell,iperm);
  t_connectiv=clock();
  cout << " : " << (t_connectiv-t_compute_graph)/(double) CLOCKS_PER_SEC << "s" << endl;
  cout.flush();

  // Familles
  cout << "Computing families    ";
  cout.flush();
  changeFamily(&meshRenum,Type,perm);
  int drv3=meshRenum.addDriver(MED_DRIVER,filename_out,meshRenum.getName());
  meshRenum.write(drv3);
  t_family=clock();
  cout << " : " << (t_family-t_connectiv)/(double) CLOCKS_PER_SEC << "s" << endl;
  cout.flush();

  // Fields
  cout << "Computing fields      ";
  cout.flush();
  bool exist_type;
  for(int ifield=0;ifield<nb_fields_tot;++ifield)
    {
      exist_type=false;
      FIELD<double> myField(MED_DRIVER,filename_in,field_names[ifield],iternumber[ifield],ordernumber[ifield]);
      FIELD<double> newField(myField);
      const SUPPORT* mySupport=newField.getSupport();
      const medGeometryElement *typesOfSupport = mySupport->getTypes();
      for(int t=0;t<mySupport->getNumberOfTypes();++t)
        {
          if(typesOfSupport[t]==Type)
            {
              exist_type=true;
              break;
            }
        }
      if(exist_type)
        {
          for(int i=0;i<mySupport->getNumberOfElements(Type);++i)
            {
              for(int j=0;j<newField.getNumberOfComponents();++j)
                {
                  newField.setValueIJ(i+1,j+1,myField.getValueIJ(iperm[i],j+1));
                }
            }
        }
      int drv=newField.addDriver(MED_DRIVER,filename_out,field_names[ifield]);
      newField.write(drv);
    }
  t_field=clock();
  cout << " : " << (t_field-t_family)/(double) CLOCKS_PER_SEC << "s" << endl;
  cout.flush();

  delete[] graph_index;
  delete[] graph;

  return 0;
}
