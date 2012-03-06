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

#include "MEDPARTITIONER_Utils.hxx"

#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"
#include "CellModel.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#ifdef HAVE_MPI2
#include <mpi.h>
#endif

using namespace std;
using namespace ParaMEDMEM;
using namespace MEDPARTITIONER;

int MEDPARTITIONER::MyGlobals::_Verbose=0;
int MEDPARTITIONER::MyGlobals::_Is0verbose=0;
int MEDPARTITIONER::MyGlobals::_Rank=(-1);
int MEDPARTITIONER::MyGlobals::_World_Size=(-1);
int MEDPARTITIONER::MyGlobals::_Randomize=0;
int MEDPARTITIONER::MyGlobals::_Atomize=0;
int MEDPARTITIONER::MyGlobals::_Creates_Boundary_Faces=0;
vector<string> MEDPARTITIONER::MyGlobals::_File_Names;
vector<string> MEDPARTITIONER::MyGlobals::_Mesh_Names;
vector<string> MEDPARTITIONER::MyGlobals::_Field_Descriptions;
vector<string> MEDPARTITIONER::MyGlobals::_General_Informations;

std::string MEDPARTITIONER::Trim(const std::string& s,const std::string& drop)
{
  string r=s;
  r.erase(r.find_last_not_of(drop)+1);
  return r.erase(0,r.find_first_not_of(drop));
}

std::string MEDPARTITIONER::IntToStr(const int i)
{
  ostringstream oss;
  oss<<i;
  return oss.str();
}

std::string MEDPARTITIONER::DoubleToStr(const double i)
{
  ostringstream oss;
  oss<<i;
  return oss.str();
}

int MEDPARTITIONER::StrToInt(const std::string& s)
{
  int res;
  istringstream iss(s);
  iss>>res;
  return res;
}

double MEDPARTITIONER::StrToDouble(const std::string& s)
{
  double res;
  istringstream iss(s);
  iss>>res;
  return res;
}

bool MEDPARTITIONER::TestArg(const char *arg, const char *argExpected, string& argValue)
{
  argValue="";
  int i;
  for (i=0; i<strlen(arg); i++)
    {
      if (arg[i]=='=') break;
      if (arg[i]!=argExpected[i]) return false;
    }
  for (int j=i+1; j<strlen(arg); j++) argValue+=arg[j];
  //cout<<"found at "<<i<<" "<<argValue<<endl;
  return true;
}

std::vector<int> MEDPARTITIONER::CreateRandomSize(const int size)
{
  vector<int> res(size);
  for (int i=0; i<size; i++) res[i]=i;
  //cvw TODO or not? srand( (unsigned)time( NULL ) );
  srand( MyGlobals::_Randomize );
  for (int i=0; i<size; i++)
    {
      int ii=rand()%size;
      int tmp=res[ii];
      res[ii]=res[i];
      res[i]=tmp;
    }
  //cout<<"CreateRandomSize "<<size<<endl;
  if (size<50) { for (int i=0; i<size; i++) cout<<res[i]<<" "; cout<<endl; }
  return res;
}

void MEDPARTITIONER::RandomizeAdj(int* xadj, int* adjncy, std::vector<int>& ran, vector<int>& vx, vector<int>& va)
//randomize a xadj and adjncy, renumbering vertices belong rand.
//work only on one processor!!!!
{
  if (MyGlobals::_World_Size>1)
    {
      cerr<<"MEDPARTITIONER::RandomizeAdj only works on one proc!"<<endl;
      return;
    }
  int size=ran.size();
  vector<int> invran(size);
  for (int i=0; i<size; i++) invran[ran[i]]=i;
  vx.resize(size+1);
  int lga=xadj[size];
  va.resize(lga);
  int jj=0;
  vx[0]=0;
  for (int i=0; i<size; i++)
    {
      int ir=ran[i];
      int ii=xadj[ir];
      int lgj=xadj[ir+1]-ii;
      for (int j=0; j<lgj; j++)
        {
          //va[jj]=adjncy[ii]; //for first debug only
          va[jj]=invran[adjncy[ii]];
          jj=jj+1;
          ii=ii+1;
        }
      vx[i+1]=jj;
    }
}

void MEDPARTITIONER::TestRandomize()
{
  //int xadj[6]={0,2,5,9,12,13}; //for first debug only
  //int adjncy[13]={1,4,0,2,4,1,3,4,2,4,4,3,4};
  int xadj[6]={0,2,5,9,12,13};
  int adjncy[13]={0,0,1,1,1,2,2,2,2,3,3,3,4};
  cout<<"TestRandomize"<<endl;
  for (int i=0; i<6; i++) cout<<xadj[i]<<" "; cout<<endl;
  for (int i=0; i<13; i++) cout<<adjncy[i]<<" "; cout<<endl;
  int size=5;
  vector<int> r=CreateRandomSize(size);
  vector<int> vx,va;
  RandomizeAdj(&xadj[0],&adjncy[0],r,vx,va);
  for (int i=0; i<vx.size(); i++) cout<<vx[i]<<" "; cout<<endl;
  for (int i=0; i<va.size(); i++) cout<<va[i]<<" "; cout<<endl;
}

std::string MEDPARTITIONER::ReprVectorOfString(const std::vector<string>& vec)
{
  if (vec.size()==0) return string(" NONE\n");
  ostringstream oss;
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i) 
    oss<<" -> '"<<*i<<"'"<<endl;
  return oss.str();
}

std::string MEDPARTITIONER::ReprVectorOfString(const std::vector<std::string>& vec, const std::string separator)
{
  if (vec.size()==0) return string(" NONE\n");
  ostringstream oss;
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i) 
    oss<<separator<<*i;
  return oss.str();
}

std::string MEDPARTITIONER::ReprMapOfStringInt(const std::map<std::string,int>& mymap)
{
  if (mymap.size()==0) return string(" NONE\n");
  ostringstream oss;
  for (map<string,int>::const_iterator i=mymap.begin(); i!=mymap.end(); ++i) 
    oss<<" -> ["<<(*i).first<<"]="<<(*i).second<<endl;
  return oss.str();
}

std::string MEDPARTITIONER::ReprMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap)
{
  if (mymap.size()==0) return string(" NONE\n");
  ostringstream oss;
  for (map< string,vector<string> >::const_iterator i=mymap.begin(); i!=mymap.end(); ++i) 
    oss<<" -> ["<<(*i).first<<"]="<<endl<<ReprVectorOfString((*i).second)<<endl;
  return oss.str();
}

std::string MEDPARTITIONER::ReprFieldDescriptions(const std::vector<std::string>& vec, const std::string separator)
{
  if (vec.size()==0) return string(" NONE\n");
  ostringstream oss;
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i)
    {
      oss<<" ->"; 
      oss<<ReprVectorOfString(DeserializeToVectorOfString(*i), separator)<<endl;
    }
  return oss.str();
}

std::string MEDPARTITIONER::SerializeFromString(const std::string& s)
//a string "hello" gives a string "    5/hello/"
//serialized_FromVectorOfString_string+SerializeFromString("toto") is
//equivalent to vector<string>.push_back("toto") on serialized_FromVectorOfString_string
{
  ostringstream oss;
  oss<<setw(5)<<s.size()<<"/"<<s<<"/";
  return oss.str();
}

std::string MEDPARTITIONER::SerializeFromVectorOfString(const std::vector<std::string>& vec)
//a vector of string gives a string
{
  ostringstream oss;
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i)
    oss<<setw(5)<<(*i).size()<<"/"<<*i<<"/";
  //string res(oss.str());
  return oss.str();
}

std::vector<std::string> MEDPARTITIONER::DeserializeToVectorOfString(const std::string& str)
//a string gives a vector of string
{
  //vector<string> res=new vector<string>;
  vector<string> res;
  size_t pos=0;
  size_t posmax=str.size();
  if (posmax==0) return res;  //empty vector
  size_t length;
  while (pos < posmax-6)  //setw(5)+" "
    {
      istringstream iss(str.substr(pos,5));
      iss>>length;
      //cout<<"length "<<length<<endl;
      if ((str[pos+5]!='/') || (str[pos+6+length]!='/'))
        {
          cerr<<"Error on string '"<<str<<"'"<<endl;;
          throw INTERP_KERNEL::Exception(LOCALIZED("Error on string"));
        }
      res.push_back(str.substr(pos+6,length));
      pos=pos+6+length+1;
    }
  return res;
}

std::string MEDPARTITIONER::EraseTagSerialized(const std::string& fromStr, const std::string& tag)
{
  vector<string> vec=DeserializeToVectorOfString(fromStr);
  vector<string> res;
  for (int i=0; i<vec.size(); i++)
    {
      if (vec[i].find(tag)==string::npos) res.push_back(vec[i]);
    }
  return MEDPARTITIONER::SerializeFromVectorOfString(res);
}

std::vector<std::string> MEDPARTITIONER::VectorizeFromMapOfStringInt(const std::map<std::string,int>& mymap)
//elements first and second of map give one elements in result vector of string
//converting formatted the int second as firsts characters ending at first slash
{
  vector<string> res;
  for (map<string,int>::const_iterator i=mymap.begin(); i!=mymap.end(); ++i)
    {
      ostringstream oss;
      oss<<(*i).second<<"/"<<(*i).first;
      res.push_back(oss.str());
    }
  return res;
}

std::map<std::string,int> MEDPARTITIONER::DevectorizeToMapOfStringInt(const std::vector<std::string>& vec)
//if existing identicals (first,second) in vector no problem, else Exception
{
  map<string,int> res;
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i)
    {
      size_t pos=0;
      size_t posmax=(*i).size();
      size_t found=(*i).find('/'); //first slash
      if ((found==string::npos) || (found<1))
        throw INTERP_KERNEL::Exception(LOCALIZED("Error aIntNumber/anyString is expected"));
      int second;
      istringstream iss((*i).substr(pos,found));
      iss>>second;
      string first=(*i).substr(pos+found+1,posmax-found);
      map<string,int>::iterator it=res.find(first);
      if (it!=res.end())
        if ((*it).second!=second)
          throw INTERP_KERNEL::Exception(LOCALIZED("Error not the same map value"));
      res[first]=second;
    }
  return res;
}

std::vector<std::string> MEDPARTITIONER::VectorizeFromMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap)
//elements first and second of map give one elements in result vector of string
//adding key map and length of second vector as first string in each serialized vector
//one serialized vector per key map
{
  vector<string> res;
  for (map< string,vector<string> >::const_iterator i=mymap.begin(); i!=mymap.end(); ++i)
    {
      vector<string> vs=(*i).second;  //a vector of string;
      ostringstream oss;
      oss<<"Keymap/"<<(*i).first<<"/"<<(*i).second.size();
      vs.insert(vs.begin(), oss.str());
      res.push_back(SerializeFromVectorOfString(vs));
    }
  return res;
}

std::map< std::string,std::vector<std::string> > MEDPARTITIONER::DevectorizeToMapOfStringVectorOfString(const std::vector<std::string>& vec)
//if existing identicals keymap in vector no problem
//duplicates in second vector
{
  map< string,vector<string> > res;
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i)
    {
      vector<string> vs=DeserializeToVectorOfString(*i);
    
      string enTete=vs[0];
      size_t posmax=enTete.size();
      size_t foundKey=enTete.find("Keymap/");
      size_t foundSizeVector=enTete.find_last_of('/');
      if ((foundKey==string::npos) || (foundKey!=0) || ((foundKey+7)>=foundSizeVector))
        throw INTERP_KERNEL::Exception(LOCALIZED("Error Keymap/anyString/aIntNumber is expected"));
      int sizeVector;
      istringstream iss(enTete.substr(foundSizeVector+1,posmax-foundSizeVector));
      iss>>sizeVector;
      string keymap=enTete.substr(foundKey+7,foundSizeVector-foundKey-7);
      //cout<<keymap<<" : sizeVector="<<enTete.substr(foundSizeVector+1,posmax-foundSizeVector)<<endl;
      for (int i=1; i<=sizeVector; i++)
        res[keymap].push_back(vs[i]); //add unconditionnaly,so merge duplicates in second vector
    }
  return res;
}

std::vector<std::string> MEDPARTITIONER::SelectTagsInVectorOfString(const std::vector<std::string>& vec, const std::string tag)
{
  vector<string> res;
  if (vec.size()==0) return res;
  //shit for unique and unique_copy for the duplicate CONSECUTIVE elements
  //I do not want to sort
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i)
    {
      if ((*i).find(tag)!=string::npos) res.push_back(*i);
    }
  return res;
}

std::vector<std::string> MEDPARTITIONER::DeleteDuplicatesInVectorOfString(const std::vector<std::string>& vec)
{
  vector<string> res;
  if (vec.size()==0) return res;
  //shit for unique and unique_copy for the duplicate CONSECUTIVE elements
  //I do not want to sort
  for (vector<string>::const_iterator i=vec.begin(); i!=vec.end(); ++i)
    {
      bool found=false;
      for (vector<string>::const_iterator j=res.begin(); j!=res.end(); ++j)
        {
          if ((*i).compare(*j)==0)
            {
              found=true;
              break;
            }
        }
      if (!found) res.push_back(*i);
    }
  return res;
}

std::map< std::string,std::vector<std::string> > MEDPARTITIONER::DeleteDuplicatesInMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap)
{
  map< string,vector<string> > res;
  for (map< string,vector<string> >::const_iterator i=mymap.begin(); i!=mymap.end(); ++i)
    res[(*i).first]=DeleteDuplicatesInVectorOfString((*i).second);
  return res;
}

/*
  void MEDPARTITIONER::SendVectorOfString(const std::vector<std::string>& vec, const int target)
  //non conseille, interblocages, utiliser sendAndReceive
  {
  throw INTERP_KERNEL::Exception(LOCALIZED("use SendAndReceiveVectorOfString please."));
  string str=SerializeFromVectorOfString(vec);
  int tag=111000;
  int size=str.length();
  MPI_Send( &size, 1, MPI_INT, target, tag, MPI_COMM_WORLD );
  MPI_Send( (void*)str.data(), str.length(), MPI_CHAR, target, tag+100, MPI_COMM_WORLD );
  cout<<"proc "<<MyGlobals::_Rank<<" : send to "<<target<<" '"<<str<<"'"<<endl;
  }
*/

/*
  std::vector<std::string> MEDPARTITIONER::RecvVectorOfString(const int source)
  //non conseille, interblocages, utiliser sendAndReceive
  {
  throw INTERP_KERNEL::Exception(LOCALIZED("use SendAndReceiveVectorOfString please."));
  int recSize=0;
  int tag=111000;
  MPI_Status status;
  MPI_Recv(&recSize, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  string recData(recSize,'x');
  MPI_Recv((void*)recData.data(), recSize, MPI_CHAR, 1, tag+100, MPI_COMM_WORLD, &status);
  //cout<<"proc "<<MyGlobals::_Rank<<" : receive from "<<source<<" '"<<recData<<"'"<<endl;
  return DeserializeToVectorOfString(recData);
  }
*/

std::vector<std::string> MEDPARTITIONER::SendAndReceiveVectorOfString(const std::vector<std::string>& vec, const int source, const int target)
//not optimized but suffisant
//return empty vector if i am not target
{
  int rank=MyGlobals::_Rank;
  
  /*for test
    ostringstream oss;
    oss<<"sendAndReceive from "<<setw(3)<<source<<" to "<<setw(3)<<target<<"-";
    string str(oss.str());
   */

  MPI_Status status;
  int tag = 111001;
  if (rank == source)
    {
      string str=SerializeFromVectorOfString(vec);
      int size=str.length();
      MPI_Send( &size, 1, MPI_INT, target, tag, MPI_COMM_WORLD );
      //cout<<"proc "<<source<<" : send "<<size<<endl;
      MPI_Send( (void*)str.data(), str.length(), MPI_CHAR, target, tag+100, MPI_COMM_WORLD );
      //cout<<"proc "<<source<<" : send    '"<<str<<"'"<<endl;
    }
  
  int recSize=0;
  if (rank == target)
    {
      MPI_Recv(&recSize, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
      //cout<<"proc "<<target<<" : receive "<<recSize<<endl;
      string recData(recSize,'x');
      MPI_Recv((void*)recData.data(), recSize, MPI_CHAR, source, tag+100, MPI_COMM_WORLD, &status);
      //cout<<"proc "<<target<<" : receive '"<<recData<<"'"<<endl;
      return DeserializeToVectorOfString(recData); //not empty one for target proc
    }
  vector<string> res;
  return res; //empty one for other proc
}

std::vector<std::string> MEDPARTITIONER::AllgathervVectorOfString(const std::vector<std::string>& vec)
//strings NO need all same size!!!!
{
  int world_size=MyGlobals::_World_Size;
  string str=SerializeFromVectorOfString(vec);
  
  /*for test
    int rank=MyGlobals::_Rank;
    ostringstream oss;
    oss<<"allgatherv from "<<setw(3)<<rank<<" to all"<<"-";
    string str(oss.str());*/
  
  vector<int> indexes(world_size);
  int size=str.length();
  MPI_Allgather(&size, 1, MPI_INT, 
                &indexes[0], 1, MPI_INT, MPI_COMM_WORLD);
  
  /*{
    ostringstream oss;
    for (int i=0; i<world_size; i++) oss<<" "<<indexes[i];
    cout<<"proc "<<rank<<" : receive '"<<oss.str()<<"'"<<endl;
    }*/
  
  //calcul of displacement
  vector< int > disp(1,0);
  for (int i=0; i<world_size; i++) disp.push_back( disp.back() + indexes[i] );
  
  string recData(disp.back(),'x');
  MPI_Allgatherv((void*)str.data(), str.length(), MPI_CHAR,
                 (void*)recData.data(), &indexes[0], &disp[0], MPI_CHAR,
                 MPI_COMM_WORLD);
  
  //really extraordinary verbose for debug
  vector<string> deserial=DeserializeToVectorOfString(recData);
  if (MyGlobals::_Verbose>1000) 
    {
      cout<<"proc "<<MyGlobals::_Rank<<" : receive '"<<recData<<"'"<<endl;
      cout<<"deserialize is : a vector of size "<<deserial.size()<<endl;
      cout<<ReprVectorOfString(deserial)<<endl;
    }
  return deserial;
}

//void MEDPARTITIONER::sendRecvVectorOfString(const std::vector<string>& vec, const int source, const int target)
//TODO

std::string MEDPARTITIONER::Cle1ToStr(const std::string& s, const int inew)
{
  ostringstream oss;
  oss<<s<<" "<<inew;
  return oss.str();
}

void MEDPARTITIONER::Cle1ToData(const std::string& cle, std::string& s, int& inew)
{
  size_t posmax=cle.size();
  size_t found=cle.find(' ');
  if ((found==string::npos) || (found<1))
    throw INTERP_KERNEL::Exception(LOCALIZED("Error 'aStringWithoutWhitespace aInt' is expected"));
  s=cle.substr(0,found);
  istringstream iss(cle.substr(found,posmax-found));
  iss>>inew;
}

std::string MEDPARTITIONER::Cle2ToStr(const std::string& s, const int inew, const int iold)
{
  ostringstream oss;
  oss<<s<<" "<<inew<<" "<<iold;
  return oss.str();
}

void MEDPARTITIONER::Cle2ToData(const std::string& cle, std::string& s, int& inew, int& iold)
{
  size_t posmax=cle.size();
  size_t found=cle.find(' ');
  if ((found==string::npos) || (found<1))
    throw INTERP_KERNEL::Exception(LOCALIZED("Error 'aStringWithoutWhitespace aInt aInt' is expected"));
  s=cle.substr(0,found);
  istringstream iss(cle.substr(found,posmax-found));
  iss>>inew>>iold;
}

std::string MEDPARTITIONER::ExtractFromDescription(const std::string& description,const std::string& tag)
{
  size_t found=description.find(tag);
  if ((found==string::npos) || (found<1))
    {
      cerr<<"ERROR : not found '"<<tag<<"' in '"<<description<<"'\n";
      throw INTERP_KERNEL::Exception(LOCALIZED("Error"));
    }
  size_t beg=found;
  size_t end=beg;
  if (description[found-1]!='/')
    {
      //find without '/'... and pray looking for first whitespace
      //something like 'idomain=0 fileName=tmp.med meshName=...'
      end=description.size();
      beg+=tag.length();
      string res=description.substr(beg,end-beg);
      found=res.find(' ');
      if (found==string::npos) found=res.length();
      res=res.substr(0,found);
      //cout<<"find without '/' !"<<tag<<"!"<<res<<"!"<<endl;
      return res;
    }
  size_t lg=StrToInt(description.substr(found-6,found));
  beg+=tag.length();
  //cout<<"find with '/' !"<<tag<<"!"<<description.substr(beg,lg-tag.length())<<"!"<<lg<<endl;
  return description.substr(beg,lg-tag.length());
}

void MEDPARTITIONER::FieldDescriptionToData(const std::string& description, 
                                            int& idomain, std::string& fileName, std::string& meshName, std::string& fieldName, int& typeField, int& DT, int& IT)
{
  idomain=StrToInt(ExtractFromDescription(description,"idomain="));
  fileName=ExtractFromDescription(description,"fileName=");
  meshName=ExtractFromDescription(description,"meshName=");
  fieldName=ExtractFromDescription(description,"fieldName=");
  typeField=StrToInt(ExtractFromDescription(description,"typeField="));
  DT=StrToInt(ExtractFromDescription(description,"DT="));
  IT=StrToInt(ExtractFromDescription(description,"IT="));
  cout<<"idomain="<<idomain<<" fileName="<<fileName<<" meshName="<<meshName<<" fieldName="<<fieldName<<endl;
}

void MEDPARTITIONER::FieldShortDescriptionToData(const std::string& description, 
                                                 std::string& fieldName, int& typeField, int& entity, int& DT, int& IT)
{
  //*meshName=ExtractFromDescription(description,"meshName=");
  fieldName=ExtractFromDescription(description,"fieldName=");
  typeField=StrToInt(ExtractFromDescription(description,"typeField="));
  entity=StrToInt(ExtractFromDescription(description,"entity="));
  DT=StrToInt(ExtractFromDescription(description,"DT="));
  IT=StrToInt(ExtractFromDescription(description,"IT="));
  //cout<<" meshName="<<*meshName<<" fieldName="<<*fieldName<<endl;
}

ParaMEDMEM::DataArrayInt* MEDPARTITIONER::CreateDataArrayIntFromVector(const std::vector<int>& v)
{
  ParaMEDMEM::DataArrayInt* p=DataArrayInt::New();
  p->alloc(v.size(),1);
  std::copy(v.begin(),v.end(),p->getPointer());
  return p;
}

ParaMEDMEM::DataArrayInt* MEDPARTITIONER::CreateDataArrayIntFromVector(const std::vector<int>& v,const int nbComponents)
{
  ParaMEDMEM::DataArrayInt* p=DataArrayInt::New();
  if (v.size()%nbComponents!=0)
    throw INTERP_KERNEL::Exception(LOCALIZED("Problem size modulo nbComponents != 0"));
  p->alloc(v.size()/nbComponents,nbComponents);
  std::copy(v.begin(),v.end(),p->getPointer());
  return p;
}

ParaMEDMEM::DataArrayDouble* MEDPARTITIONER::CreateDataArrayDoubleFromVector(const std::vector<double>& v)
{
  ParaMEDMEM::DataArrayDouble* p=DataArrayDouble::New();
  p->alloc(v.size(),1);
  std::copy(v.begin(),v.end(),p->getPointer());
  return p;
}

std::vector<std::string> MEDPARTITIONER::BrowseFieldDouble(const MEDCouplingFieldDouble* fd)
//quick almost human readable information on a field double

/*
  example done by fd->simpleRepr() :
  FieldDouble with name : "VectorFieldOnCells"
  Description of field is : ""
  FieldDouble space discretization is : P0
  FieldDouble time discretization is : One time label. Time is defined by iteration=0 order=1 and time=2.
  Time unit is : ""
  FieldDouble nature of field is : NoNature
  FieldDouble default array has 3 components and 30000 tuples.
  FieldDouble default array has following info on components : "vx" "vy" "vz"
  Mesh support information :
  __________________________
  Unstructured mesh with name : "testMesh"
  Description of mesh : ""
  Time attached to the mesh [unit] : 0 []
  Iteration : -1 Order : -1
  Mesh dimension : 3
  Space dimension : 3
  Info attached on space dimension : "" "" ""
  Number of nodes : 33201
  Number of cells : 30000
  Cell types present : NORM_HEXA8
*/

{
  vector<string> res;
  //res.push_back("fieldName="); res.back()+=fd->getName();
  //not saved in file? res.push_back("fieldDescription="); res.back()+=fd->getDescription();
  //ret << "FieldDouble space discretization is : " << _type->getStringRepr() << "\n";
  //ret << "FieldDouble time discretization is : " << _time_discr->getStringRepr() << "\n";
  //ret << "FieldDouble nature of field is : " << MEDCouplingNatureOfField::getRepr(_nature) << "\n";
  if (fd->getArray())
    {
      int nb=fd->getArray()->getNumberOfComponents();
      res.push_back("nbComponents="); res.back()+=IntToStr(nb);
      //ret << "FieldDouble default array has " << nbOfCompo << " components and " << getArray()->getNumberOfTuples() << " tuples.\n";
      //ret << "FieldDouble default array has following info on components : ";
      for (int i=0; i<nb; i++)
        {
          //ret << "\"" << getArray()->getInfoOnComponent(i) << "\" ";
          res.push_back("componentInfo");
          res.back()+=IntToStr(i)+"="+fd->getArray()->getInfoOnComponent(i);
        }
    }
  else
    {
      res.push_back("nbComponents=0");  //unknown
    }
  return res;
}

std::vector<std::string> MEDPARTITIONER::BrowseAllFields(const std::string& myfile)
//quick almost human readable information on all fields in a .med file
{
  vector<string> res;
  vector<string> meshNames=MEDLoader::GetMeshNames(myfile.c_str());
  
  for (int i=0; i<meshNames.size(); i++)
    {
      vector<string> fieldNames=
        MEDLoader::GetAllFieldNamesOnMesh(myfile.c_str(),meshNames[i].c_str());
      for (int j = 0; j < fieldNames.size(); j++)
        {
          vector< ParaMEDMEM::TypeOfField > typeFields=
            MEDLoader::GetTypesOfField(myfile.c_str(), meshNames[i].c_str(), fieldNames[j].c_str());
          for (int k = 0; k < typeFields.size(); k++)
            {
              vector< pair< int, int > > its=
                MEDLoader::GetFieldIterations(typeFields[k], myfile.c_str(), meshNames[i].c_str(), fieldNames[j].c_str());
              if (MyGlobals::_Is0verbose>100) cout<<"fieldName "<<fieldNames[j].c_str()<<" typeField "<<typeFields[k]<<" its.size() "<<its.size()<<endl;
              for (int m = 0; m < its.size(); m++)
                {
                  vector<string> resi;
                  resi.push_back("fileName="); resi.back()+=myfile;
                  resi.push_back("meshName="); resi.back()+=meshNames[i];
                  resi.push_back("fieldName="); resi.back()+=fieldNames[j];
                  resi.push_back("typeField="); resi.back()+=IntToStr((int)typeFields[k]);
                  resi.push_back("DT="); resi.back()+=IntToStr((int)its[m].first);
                  resi.push_back("IT="); resi.back()+=IntToStr((int)its[m].second);
                  res.push_back(SerializeFromVectorOfString(resi));
                }
            }
        }
    }
  return res;
}

std::vector<std::string> MEDPARTITIONER::GetInfosOfField(const char *fileName, const char *meshName, const int idomain)
{
  const int lggeom=10;
  const med_geometry_type GEOMTYPE[lggeom]={ //MED_N_CELL_FIXED_GEO] = { 
    //MED_POINT1,
    //MED_SEG2,
    //MED_SEG3,
    //MED_SEG4,
    //MED_TRIA3,
    //MED_QUAD4,
    //MED_TRIA6,
    //MED_TRIA7,
    //MED_QUAD8,
    //MED_QUAD9,
    MED_TETRA4,
    MED_PYRA5,
    MED_PENTA6,
    MED_HEXA8,
    MED_OCTA12,
    MED_TETRA10,
    MED_PYRA13,
    MED_PENTA15,
    MED_HEXA20,
    MED_HEXA27,
    //MED_POLYGON,
    //MED_POLYHEDRON 
  };

  const char * const GEOMTYPENAME[lggeom]={
    //"MED_POINT1",
    //"MED_SEG2",
    //"MED_SEG3",
    //"MED_SEG4",
    //"MED_TRIA3",
    //"MED_QUAD4",
    //"MED_TRIA6",
    //"MED_TRIA7",
    //"MED_QUAD8",
    //"MED_QUAD9",
    "MED_TETRA4",
    "MED_PYRA5",
    "MED_PENTA6",
    "MED_HEXA8",
    "MED_OCTA12",
    "MED_TETRA10",
    "MED_PYRA13",
    "MED_PENTA15",
    "MED_HEXA20",
    "MED_HEXA27",
    //"MED_POLYGONE",
    //"MED_POLYEDRE",
  };


  const int lgentity=3;
  const med_entity_type ENTITYTYPE[lgentity]={ //MED_N_ENTITY_TYPES+2]={
    //MED_UNDEF_ENTITY_TYPE,
    MED_CELL,
    //MED_DESCENDING_FACE,
    //MED_DESCENDING_EDGE,
    MED_NODE,
    MED_NODE_ELEMENT,
    //MED_STRUCT_ELEMENT,
    //MED_UNDEF_ENTITY_TYPE
  };

  const char * const ENTITYTYPENAME[lgentity]={ //MED_N_ENTITY_TYPES+2]={
    //"MED_UNDEF_ENTITY_TYPE",
    "MED_CELL",
    //"MED_FACE",
    //"MED_ARETE",
    "MED_NODE",
    "MED_NODE_ELEMENT",
    //"MED_STRUCT_ELEMENT",
    //"MED_UNDEF_ENTITY_TYPE"
  };
  
  vector<string> res;
  med_idt fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  if (MyGlobals::_Verbose>20) cout<<"on filename "<<fileName<<" nbOfField "<<nbFields<<endl;
  //
  med_field_type typcha;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  for(int i=1; i<=nbFields; i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> dt_unit=new char[MED_LNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      string curMeshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE+1);
      for (int k=1; k<=nbPdt; k++)
        {
          MEDfieldComputingStepInfo(fid,nomcha,k,&numdt,&numo,&dt);
          if (MyGlobals::_Verbose>20) 
            cout<<"on filename "<<fileName<<" field "<<i<<" fieldName "<<curFieldName<<" meshName "<<curMeshName<<
              " typ "<<typcha<<" nbComponent "<<ncomp<<" nbPdt "<<nbPdt<<" noPdt "<<k<<
              " ndt "<<numdt<<" nor "<<numo<<" dt "<<dt<<endl;
          for (int ie=0; ie<lgentity; ie++)
            {
              for (int j=0; j<lggeom; j++)
                {
                  int profilesize=0,nbi=0;
                  med_entity_type enttype=ENTITYTYPE[ie];
                  //enttype=MED_NODE;enttype=MED_CELL;enttype=MED_NODE_ELEMENT;
                  char pflname[MED_NAME_SIZE+1]="";
                  char locname[MED_NAME_SIZE+1]="";
                  med_int nbofprofile=MEDfieldnProfile(fid,nomcha,numdt,numo,enttype,GEOMTYPE[j],pflname,locname);

                  /*
                    med_proto.h firefox file:///export/home/.../med-3.0.3/doc/html.dox/index.html
                    MEDfieldnValueWithProfileByName(const med_idt fid, const char * const fieldname,const med_int numdt,const med_int numit,
                    const med_entity_type entitype, const med_geometry_type geotype, const char * const profilename,
                    const med_storage_mode storagemode,med_int * const profilesize,
                    char * const localizationname, med_int * const nbofintegrationpoint);

                    MEDfieldnValueWithProfile(const med_idt fid, const char * const fieldname,const med_int numdt,const med_int numit,
                    const med_entity_type entitype, const med_geometry_type geotype,
                    const int profileit,
                    const med_storage_mode storagemode,char * const profilename ,med_int * const profilesize,
                    char * const localizationname, med_int * const nbofintegrationpoint);
                  */
                  int profileit=1;
                  if (enttype==MED_NODE)
                    {
                      med_geometry_type mygeomtype=MED_UNDEF_ENTITY_TYPE;
                      med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,enttype,mygeomtype,profileit,
                                                                MED_COMPACT_PFLMODE,pflname,&profilesize,locname,&nbi);
                      if (nbOfVal>0)
                        {
                          if (MyGlobals::_Verbose>20)
                            cout<<"on filename "<<fileName<<" entity "<<enttype<<" nbOfVal with "<<
                              nbofprofile<<" profile(s) for geomType (AUCUN) nbOfVal "<<
                              nbOfVal<<" profilName '"<<pflname<<"' profileSize "<<profilesize<<" nbPtGauss "<<nbi<<endl;
                          vector<string> resi;
                          resi.push_back("idomain="); resi.back()+=IntToStr(idomain);
                          resi.push_back("fileName="); resi.back()+=fileName;
                          resi.push_back("meshName="); resi.back()+=curMeshName;
                          resi.push_back("fieldName="); resi.back()+=curFieldName;
                          resi.push_back("typeField="); resi.back()+=IntToStr((int)ON_NODES);
                          resi.push_back("typeData="); resi.back()+=IntToStr((int)typcha);  //6 for double?
                          resi.push_back("nbComponent="); resi.back()+=IntToStr((int)ncomp);
                          resi.push_back("DT="); resi.back()+=IntToStr((int)numdt);
                          resi.push_back("IT="); resi.back()+=IntToStr((int)numo);
                          resi.push_back("time="); resi.back()+=DoubleToStr(dt);
                          resi.push_back("entity="); resi.back()+=IntToStr((int)enttype);
                          resi.push_back("entityName="); resi.back()+=ENTITYTYPENAME[ie];
                          resi.push_back("nbOfVal="); resi.back()+=IntToStr((int)nbOfVal);
                          resi.push_back("profilName="); resi.back()+=pflname;
                          resi.push_back("profileSize="); resi.back()+=IntToStr((int)profilesize);
                          resi.push_back("nbPtGauss="); resi.back()+=IntToStr((int)nbi);
                          res.push_back(SerializeFromVectorOfString(resi));
                        }
                      break; //on nodes no need to scute all geomtype
                    }
                  else
                    {
                      med_geometry_type mygeomtype=GEOMTYPE[j];
                      med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,enttype,mygeomtype,profileit,
                                                                MED_COMPACT_PFLMODE,pflname,&profilesize,locname,&nbi);
                      if (nbOfVal>0)
                        {
                          if (MyGlobals::_Verbose>20)
                            cout<<"on filename "<<fileName<<" entity "<<enttype<<" nbOfVal with "<<
                              nbofprofile<<" profile(s) for geomType "<<
                              GEOMTYPE[j]<<" "<<GEOMTYPENAME[j]<<" nbOfVal "<<
                              nbOfVal<<" profilName '"<<pflname<<"' profileSize "<<profilesize<<" nbPtGauss "<<nbi<<endl;
                          int typeField=(-1); //unknown
                          if (enttype==MED_CELL) typeField=ON_CELLS;
                          if (enttype==MED_NODE_ELEMENT) typeField=ON_GAUSS_NE;
                          //if (enttype==??) typeField=ON_GAUSS_PT;
                          vector<string> resi;
                          resi.push_back("idomain="); resi.back()+=IntToStr(idomain);
                          resi.push_back("fileName="); resi.back()+=fileName;
                          resi.push_back("meshName="); resi.back()+=curMeshName;
                          resi.push_back("fieldName="); resi.back()+=curFieldName;
                          resi.push_back("typeField="); resi.back()+=IntToStr((int)typeField);
                          resi.push_back("typeData="); resi.back()+=IntToStr((int)typcha);  //6 for double?
                          resi.push_back("nbComponent="); resi.back()+=IntToStr((int)ncomp);
                          resi.push_back("DT="); resi.back()+=IntToStr((int)numdt);
                          resi.push_back("IT="); resi.back()+=IntToStr((int)numo);
                          resi.push_back("time="); resi.back()+=DoubleToStr(dt);
                          resi.push_back("entity="); resi.back()+=IntToStr((int)enttype);
                          resi.push_back("entityName="); resi.back()+=ENTITYTYPENAME[ie];
                          resi.push_back("geomType="); resi.back()+=IntToStr((int)GEOMTYPE[j]);
                          resi.push_back("geomTypeName="); resi.back()+=GEOMTYPENAME[j];
                          resi.push_back("nbOfVal="); resi.back()+=IntToStr((int)nbOfVal);
                          resi.push_back("profilName="); resi.back()+=pflname;
                          resi.push_back("profileSize="); resi.back()+=IntToStr((int)profilesize);
                          resi.push_back("nbPtGauss="); resi.back()+=IntToStr((int)nbi);
                          if (typeField==(-1))
                            {
                              cout<<"WARNING : unknown typeField for entity type "<<enttype<<endl<<
                                SerializeFromVectorOfString(resi)<<endl;
                              continue;  //do not register push_back
                            }
                          res.push_back(SerializeFromVectorOfString(resi));
                        }
                    }
                }
            }
        }
    }
  delete [] maa_ass;
  delete [] nomcha;
  MEDfileClose(fid);
  if (MyGlobals::_Verbose>10) cout<<"detected fields:\n"<<ReprVectorOfString(res)<<endl;
  return res;
}

std::vector<std::string> MEDPARTITIONER::BrowseAllFieldsOnMesh(const std::string& myfile, const std::string& mymesh, const int idomain)
//quick almost human readable information on all fields on a mesh in a .med file
{
  vector<string> res=GetInfosOfField(myfile.c_str(),mymesh.c_str(),idomain);
  return res;
  
  /*obsolete do no work on GetTypesOfField ON_GAUSS_NE
    vector<string> res;
    vector<string> meshNames;
    meshNames.push_back(mymesh);
  
    for (int i=0; i<meshNames.size(); i++)
    {
    vector<string> fieldNames=
    MEDLoader::GetAllFieldNamesOnMesh(myfile.c_str(),meshNames[i].c_str());
    for (int j=0; j<fieldNames.size(); j++)
    {
    vector< ParaMEDMEM::TypeOfField > typeFields=
    MEDLoader::GetTypesOfField(myfile.c_str(), meshNames[i].c_str(), fieldNames[j].c_str());
    //if (MyGlobals::_Is0verbose>100) cout<<"fieldName "<<fieldNames[j].c_str()<<"typeField.size "<<typeFields.size()<<endl;
    if (typeFields.size()==0) {
    cerr<<"problem : fieldNames "<<fieldNames[j]<<" without GetTypesOfField ! Type of field specified not managed"<<endl;
    //typeFields.push_back(ON_GAUSS_NE);
    }
    for (int k=0; k<typeFields.size(); k++)
    {
    vector< pair< int, int > > its;
    its=MEDLoader::GetFieldIterations(typeFields[k], myfile.c_str(), meshNames[i].c_str(), fieldNames[j].c_str());
    //if (typeFields[k]==ON_GAUSS_NE) its.push_back(make_pair(5,6));
    if (MyGlobals::_Is0verbose>100) cout<<"fieldName "<<fieldNames[j].c_str()<<" typeField "<<typeFields[k]<<" its.size() "<<its.size()<<endl;
    for (int m = 0; m < its.size(); m++)
    {
    vector<string> resi;
    resi.push_back("idomain="); resi.back()+=IntToStr(idomain);
    resi.push_back("fileName="); resi.back()+=myfile;
    resi.push_back("meshName="); resi.back()+=meshNames[i];
    resi.push_back("fieldName="); resi.back()+=fieldNames[j];
    resi.push_back("typeField="); resi.back()+=IntToStr((int)typeFields[k]);
    resi.push_back("DT="); resi.back()+=IntToStr((int)its[m].first);
    resi.push_back("IT="); resi.back()+=IntToStr((int)its[m].second);
    //cout<<"BrowseAllFieldsOnMesh add "<<resi.size()<<endl;
    res.push_back(SerializeFromVectorOfString(resi));
    }
    }
    }
    }
    return res;
  */
}

/*!
  Sends content of \a vec to processor \a target. To be used with \a RecvDoubleVec method.
  \param vec vector to be sent
  \param target processor id of the target
*/
void MEDPARTITIONER::SendDoubleVec(const std::vector<double>& vec, const int target)
{
  int tag = 111002;
  int size=vec.size();
  if (MyGlobals::_Verbose>1000) 
    cout<<"proc "<<MyGlobals::_Rank<<" : --> SendDoubleVec "<<size<<endl;
#ifdef HAVE_MPI2
  MPI_Send(&size, 1, MPI_INT, target, tag, MPI_COMM_WORLD);
  MPI_Send(const_cast<double*>(&vec[0]), size, MPI_DOUBLE, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill vector<int> vec.
  To be used with \a SendDoubleVec method.

  \param vec vector that is filled
  \param source processor id of the incoming messages
*/
std::vector<double>* MEDPARTITIONER::RecvDoubleVec(const int source)
{
  int tag = 111002;
  int size;
#ifdef HAVE_MPI2
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000) 
    cout<<"proc "<<MyGlobals::_Rank<<" : <-- RecvDoubleVec "<<size<<endl;;
  vector<double>* vec=new vector<double>;
  vec->resize(size);
  MPI_Recv(&vec[0], size, MPI_DOUBLE, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return vec;
}

void MEDPARTITIONER::RecvDoubleVec(std::vector<double>& vec, const int source)
{
  int tag = 111002;
  int size;
#ifdef HAVE_MPI2
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    cout<<"proc "<<MyGlobals::_Rank<<" : <-- RecvDoubleVec "<<size<<endl;;
  vec.resize(size);
  MPI_Recv(&vec[0], size, MPI_DOUBLE, source, tag+100, MPI_COMM_WORLD, &status);
#endif
}
/*!
  Sends content of \a vec to processor \a target. To be used with \a RecvIntVec method.
  \param vec vector to be sent
  \param target processor id of the target
*/
void MEDPARTITIONER::SendIntVec(const std::vector<int>& vec, const int target)
{
  int tag = 111003;
  int size=vec.size();
  if (MyGlobals::_Verbose>1000)
    cout<<"proc "<<MyGlobals::_Rank<<" : --> SendIntVec "<<size<<endl; //cvw 
#ifdef HAVE_MPI2
  MPI_Send(&size, 1, MPI_INT, target, tag, MPI_COMM_WORLD);
  MPI_Send(const_cast<int*>(&vec[0]), size,MPI_INT, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill vector<int> vec.
  To be used with \a SendIntVec method.
  \param vec vector that is filled
  \param source processor id of the incoming messages
*/
std::vector<int>* MEDPARTITIONER::RecvIntVec(const int source)
{
  int tag = 111003;
  int size;
#ifdef HAVE_MPI2
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000) 
    cout<<"proc "<<MyGlobals::_Rank<<" : <-- RecvIntVec "<<size<<endl; //cvw 
  vector<int>* vec=new vector<int>;
  vec->resize(size);
  MPI_Recv(&vec[0], size, MPI_INT, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return vec;
}

void MEDPARTITIONER::RecvIntVec(std::vector<int>& vec, const int source)
{
  int tag = 111003;
  int size;
#ifdef HAVE_MPI2
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    cout<<"proc "<<MyGlobals::_Rank<<" : <-- RecvIntVec "<<size<<endl; //cvw 
  vec.resize(size);
  MPI_Recv(&vec[0], size, MPI_INT, source, tag+100, MPI_COMM_WORLD,&status);
#endif
}

/*!
  Sends content of \a dataArrayInt to processor \a target. 
  To be used with \a RecvDataArrayInt method.
  \param da dataArray to be sent
  \param target processor id of the target
*/
void MEDPARTITIONER::SendDataArrayInt(const ParaMEDMEM::DataArrayInt* da, const int target)
{
  if (da==0) throw INTERP_KERNEL::Exception(LOCALIZED("Problem send DataArrayInt* NULL"));
  int tag = 111004;
  int size[3];
  size[0]=da->getNbOfElems();
  size[1]=da->getNumberOfTuples();
  size[2]=da->getNumberOfComponents();
  if (MyGlobals::_Verbose>1000) 
    cout<<"proc "<<MyGlobals::_Rank<<" : --> SendDataArrayInt "<<size[0]<<endl; //cvw 
#ifdef HAVE_MPI2
  MPI_Send(&size, 3, MPI_INT, target, tag, MPI_COMM_WORLD);
  const int * p=da->getConstPointer();
  MPI_Send(const_cast<int*>(&p[0]), size[0] ,MPI_INT, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill dataArrayInt da.
  To be used with \a SendIntVec method.
  \param da dataArrayInt that is filled
  \param source processor id of the incoming messages
*/
ParaMEDMEM::DataArrayInt* MEDPARTITIONER::RecvDataArrayInt(const int source)
{
  int tag = 111004;
  int size[3];
#ifdef HAVE_MPI2
  MPI_Status status;
  MPI_Recv(size, 3, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    cout<<"proc "<<MyGlobals::_Rank<<" : <-- RecvDataArrayInt "<<size[0]<<endl; //cvw 
  if (size[0]!=(size[1]*size[2]))
    throw INTERP_KERNEL::Exception(LOCALIZED("Problem in RecvDataArrayInt incoherent sizes"));
  ParaMEDMEM::DataArrayInt* da=ParaMEDMEM::DataArrayInt::New();
  da->alloc(size[1],size[2]);
  int * p=da->getPointer();
  MPI_Recv(const_cast<int*>(&p[0]), size[0], MPI_INT, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return da;
}

/*!
  Sends content of \a dataArrayInt to processor \a target. 
  To be used with \a RecvDataArrayDouble method.
  \param da dataArray to be sent
  \param target processor id of the target
*/
void MEDPARTITIONER::SendDataArrayDouble(const ParaMEDMEM::DataArrayDouble* da, const int target)
{
  if (da==0) throw INTERP_KERNEL::Exception(LOCALIZED("Problem send DataArrayDouble* NULL"));
  int tag = 111005;
  int size[3];
  size[0]=da->getNbOfElems();
  size[1]=da->getNumberOfTuples();
  size[2]=da->getNumberOfComponents();
  if (MyGlobals::_Verbose>1000) 
    cout<<"proc "<<MyGlobals::_Rank<<" : --> SendDataArrayDouble "<<size[0]<<endl; //cvw 
#ifdef HAVE_MPI2
  MPI_Send(&size, 3, MPI_INT, target, tag, MPI_COMM_WORLD);
  const double * p=da->getConstPointer();
  MPI_Send(const_cast<double*>(&p[0]), size[0] ,MPI_DOUBLE, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill dataArrayDouble da.
  To be used with \a SendDoubleVec method.
  \param da dataArrayDouble that is filled
  \param source processor id of the incoming messages
*/
ParaMEDMEM::DataArrayDouble* MEDPARTITIONER::RecvDataArrayDouble(const int source)
{
  int tag = 111005;
  int size[3];
#ifdef HAVE_MPI2
  MPI_Status status;
  MPI_Recv(size, 3, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    cout<<"proc "<<MyGlobals::_Rank<<" : <-- RecvDataArrayDouble "<<size[0]<<endl; //cvw 
  if (size[0]!=(size[1]*size[2]))
    throw INTERP_KERNEL::Exception(LOCALIZED("Problem in RecvDataArrayDouble incoherent sizes"));
  ParaMEDMEM::DataArrayDouble* da=ParaMEDMEM::DataArrayDouble::New();
  da->alloc(size[1],size[2]);
  double * p=da->getPointer();
  MPI_Recv(const_cast<double*>(&p[0]), size[0], MPI_DOUBLE, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return da;
}

ParaMEDMEM::MEDCouplingUMesh* MEDPARTITIONER::CreateEmptyMEDCouplingUMesh()
//create empty MEDCouplingUMesh* dim 3
{
  ParaMEDMEM::MEDCouplingUMesh* umesh=ParaMEDMEM::MEDCouplingUMesh::New();
  umesh->setMeshDimension(3);
  umesh->allocateCells(0);
  umesh->finishInsertingCells();
  ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
  myCoords->alloc(0,3);
  umesh->setCoords(myCoords);
  umesh->setName("EMPTY");
  myCoords->decrRef();
  umesh->checkCoherency();
  return umesh;
}

void MEDPARTITIONER::TestVectorOfStringMpi()
{
  int rank=MyGlobals::_Rank;
  int world_size=MyGlobals::_World_Size;
  vector<string> myVector;
  ostringstream oss;
  oss<<"hello from "<<setw(5)<<rank<<" "<<string(rank+1,'n')<<
    " next is an empty one";
  myVector.push_back(oss.str());
  myVector.push_back("");
  myVector.push_back("next is an singleton");
  myVector.push_back("1");
  
  if (rank==0)
    {
      /*
        string s0=SerializeFromVectorOfString(myVector);
        cout<<"s0 is : a string '"<<s0<<"'"<<endl;
        vector<string> v0=DeserializeToVectorOfString(s0);
        cout<<"v0 is : a vector of size "<<v0.size()<<endl;
        cout<<ReprVectorOfString(v0)<<endl;
      */
      string s0=SerializeFromVectorOfString(myVector);
      vector<string> res=DeserializeToVectorOfString(s0);
      if (res.size()!=myVector.size()) 
        throw INTERP_KERNEL::Exception(LOCALIZED("Problem in (de)serialise VectorOfString incoherent sizes"));
      for (int i=0; i<myVector.size(); i++)
        if (res[i]!=myVector[i])
          throw INTERP_KERNEL::Exception(LOCALIZED("Problem in (de)serialise VectorOfString incoherent elements"));
    }

  for (int i=0; i<world_size; i++)
    {
      for (int j=0; j<world_size; j++)
        {
          vector<string> res=SendAndReceiveVectorOfString(myVector, i, j);
          if ((rank==j) && MyGlobals::_Verbose>20)
            cout<<"proc "<<rank<<" : receive \n"<<ReprVectorOfString(res)<<endl;
          if (rank==j)
            {
              if (res.size()!=myVector.size()) 
                throw INTERP_KERNEL::Exception(LOCALIZED("Problem in SendAndReceiveVectorOfString incoherent sizes"));
              for (int i=1; i<myVector.size(); i++) //first is different
                if (res[i]!=myVector[i])
                  throw INTERP_KERNEL::Exception(LOCALIZED("Problem in SendAndReceiveVectorOfString incoherent elements"));
            }
          else 
            {
              if (res.size()!=0) 
                throw INTERP_KERNEL::Exception(LOCALIZED("Problem in SendAndReceiveVectorOfString size have to be 0"));
            }
        }
    }
  vector<string> res=AllgathervVectorOfString(myVector);
  //sometimes for test
  res=AllgathervVectorOfString(myVector);
  res=AllgathervVectorOfString(myVector);
  if (rank==0 && MyGlobals::_Verbose>20)
    cout<<"proc "<<rank<<" : receive \n"<<ReprVectorOfString(res)<<endl;
  if (res.size()!=myVector.size()*world_size) 
    throw INTERP_KERNEL::Exception(LOCALIZED("Problem in AllgathervVectorOfString incoherent sizes"));
  int jj=(-1);
  for (int j=0; j<world_size; j++)
    {
      for (int i=0; i<myVector.size(); i++)
        {
          jj=jj+1;
          if (i==0) continue; //first is different
          if (res[jj]!=myVector[i])
            throw INTERP_KERNEL::Exception(LOCALIZED("Problem in AllgathervVectorOfString incoherent elements"));
        }
    }
  if (MyGlobals::_Verbose) cout<<"proc "<<rank<<" : OK TestVectorOfStringMpi END"<< endl;
}

void MEDPARTITIONER::TestMapOfStringIntMpi()
{
  int rank=MyGlobals::_Rank;
  //int world_size=MyGlobals::_World_Size;
  map<string,int> myMap;
  myMap["one"]=1;
  myMap["two"]=22;  //a bug
  myMap["three"]=3;
  myMap["two"]=2; //last speaking override
  
  if (rank==0)
    {
      vector<string> v2=VectorizeFromMapOfStringInt(myMap);
      /*
        cout<<"v2 is : a vector of size "<<v2.size()<<endl;
        cout<<ReprVectorOfString(v2)<<endl;
      */
      map<string,int> m3=DevectorizeToMapOfStringInt(v2);
      if (ReprMapOfStringInt(m3)!=ReprMapOfStringInt(myMap))
        throw INTERP_KERNEL::Exception(LOCALIZED("Problem in (de)vectorize MapOfStringInt"));
    }
    
  vector<string> v2=AllgathervVectorOfString(VectorizeFromMapOfStringInt(myMap));
  if (rank==0 && MyGlobals::_Verbose>20)
    {
      cout<<"v2 is : a vector of size "<<v2.size()<<endl;
      cout<<ReprVectorOfString(v2)<<endl;
      map<string,int> m2=DevectorizeToMapOfStringInt(v2);
      cout<<"m2 is : a map of size "<<m2.size()<<endl;
      cout<<ReprMapOfStringInt(m2)<<endl;
    }
  if (MyGlobals::_Verbose) cout<<"proc "<<rank<<" : OK TestMapOfStringIntMpi END"<< endl;
}

void MEDPARTITIONER::TestMapOfStringVectorOfStringMpi()
{
  int rank=MyGlobals::_Rank;
  //int world_size=MyGlobals::_World_Size;
  vector<string> myVector;
  ostringstream oss;
  oss<<"hello from "<<setw(5)<<MyGlobals::_Rank<<" "<<string(rank+1,'n')<<
    " next is an empty one";
  myVector.push_back(oss.str());
  myVector.push_back("");
  myVector.push_back("next is an singleton");
  myVector.push_back("1");
  
  if (rank==0)
    {
      map< string,vector<string> > m2;
      m2["first key"]=myVector;
      m2["second key"]=myVector;
      vector<string> v2=VectorizeFromMapOfStringVectorOfString(m2);
      map< string,vector<string> > m3=DevectorizeToMapOfStringVectorOfString(v2);
      if (rank==0 && MyGlobals::_Verbose>20)
        cout<<"m2 is : a MapOfStringVectorOfString of size "<<m2.size()<<endl;
      cout<<ReprMapOfStringVectorOfString(m2)<<endl;
      cout<<"v2 is : a vector of size "<<v2.size()<<endl;
      cout<<ReprVectorOfString(v2)<<endl;
      cout<<"m3 is : a map of size "<<m3.size()<<endl;
      cout<<ReprMapOfStringVectorOfString(m3)<<endl;
      if (ReprMapOfStringVectorOfString(m3)!=ReprMapOfStringVectorOfString(m2))
        throw INTERP_KERNEL::Exception(LOCALIZED("Problem in (de)vectorize MapOfStringVectorOfString"));
    }
    
  map< string,vector<string> > m4;
  m4["1rst key"]=myVector;
  m4["2snd key"]=myVector;
  vector<string> v4=AllgathervVectorOfString(VectorizeFromMapOfStringVectorOfString(m4));
  if (rank==0 && MyGlobals::_Verbose>20)
    {
      map< string,vector<string> > m5=DevectorizeToMapOfStringVectorOfString(v4);
      map< string,vector<string> > m6=DeleteDuplicatesInMapOfStringVectorOfString(m5);
      cout<<"m5 is : a map of size "<<m5.size()<<endl;
      cout<<ReprMapOfStringVectorOfString(m5)<<endl;
      cout<<"m6 is : a map from m5 with deleteDuplicates of size "<<m6.size()<<endl;
      cout<<ReprMapOfStringVectorOfString(m6)<<endl;
    }
  if (MyGlobals::_Verbose) cout<<"proc "<<rank<<" : OK TestMapOfStringVectorOfStringMpi END"<< endl;
}

void MEDPARTITIONER::TestDataArrayMpi()
{
  int rank=MyGlobals::_Rank;
  //int
  {
    ParaMEDMEM::DataArrayInt* send=ParaMEDMEM::DataArrayInt::New();
    ParaMEDMEM::DataArrayInt* recv=0;
    int nbOfTuples=5;
    int numberOfComponents=3;
    send->alloc(nbOfTuples,numberOfComponents);
    vector<int> vals;
    for (int j=0; j<nbOfTuples; j++)
      for (int i=0; i<numberOfComponents; i++) vals.push_back((j+1)*10+i+1);
    std::copy(vals.begin(),vals.end(),send->getPointer());
    if (rank==0) SendDataArrayInt(send, 1);
    if (rank==1) recv=RecvDataArrayInt(0);
    if (rank==1 && MyGlobals::_Verbose>20)
      {
        cout<<send->repr()<<endl;
        cout<<recv->repr()<<endl;
      }
    if (rank==1)
      {
        if (send->repr()!=recv->repr())
          throw INTERP_KERNEL::Exception(LOCALIZED("Problem in send&recv DataArrayInt"));
      }
    send->decrRef();
    if (rank==1) recv->decrRef();
  }
  //double
  {
    ParaMEDMEM::DataArrayDouble* send=ParaMEDMEM::DataArrayDouble::New();
    ParaMEDMEM::DataArrayDouble* recv=0;
    int nbOfTuples=5;
    int numberOfComponents=3;
    send->alloc(nbOfTuples,numberOfComponents);
    vector<double> vals;
    for (int j=0; j<nbOfTuples; j++)
      for (int i=0; i<numberOfComponents; i++) vals.push_back(double(j+1)+double(i+1)/10);
    std::copy(vals.begin(),vals.end(),send->getPointer());
    if (rank==0) SendDataArrayDouble(send, 1);
    if (rank==1) recv=RecvDataArrayDouble(0);
    if (rank==1 && MyGlobals::_Verbose>20)
      {
        cout<<send->repr()<<endl;
        cout<<recv->repr()<<endl;
      }
    if (rank==1)
      {
        if (send->repr()!=recv->repr())
          throw INTERP_KERNEL::Exception(LOCALIZED("Problem in send&recv DataArrayDouble"));
      }
    send->decrRef();
    if (rank==1) recv->decrRef();
  }
  
  if (MyGlobals::_Verbose) cout<<"proc "<<rank<<" : OK TestDataArrayMpi END"<< endl;
}

void MEDPARTITIONER::TestPersistantMpi0To1(int taille, int nb)
{
  double temps_debut=MPI_Wtime();
  int rang=MyGlobals::_Rank;
  vector<int> x, y;
  int tag=111111;
  MPI_Request requete0, requete1;
  MPI_Status statut;
  int ok=0;
  string res;
  if (rang==0)
    {
      x.resize(taille);
      MPI_Ssend_init(&x[0], taille, MPI_INT, 1, tag, MPI_COMM_WORLD , &requete0);
      for(int k=0; k<nb; k++)
        {
          for (int i=0; i<taille; ++i) x[i]=k;
          //Envoi d’un gros message --> cela peut prendre du temps
          MPI_Start(&requete0);
          //Traitement sequentiel independant de "x"
          //...
          MPI_Wait(&requete0, &statut);
          //Traitement sequentiel impliquant une modification de "x" en memoire
          //x=...
        }
      MPI_Request_free(&requete0);
    }
  else if (rang == 1)
    {
      y.resize(taille);
      MPI_Recv_init(&y[0], taille,  MPI_INT, 0, tag, MPI_COMM_WORLD , &requete1);
      for(int k=0; k<nb; k++)
        {
          //Pre-traitement sequentiel
          //...
          for (int i=0; i<taille; ++i) y[i]=(-1);
          //Reception du gros message --> cela peut prendre du temps
          MPI_Start(&requete1);
          //Traitement sequentiel independant de "y"
          //...
          MPI_Wait(&requete1, &statut);
          //Traitement sequentiel dependant de "y"
          //...=f(y)
          int nb=0;
          for (int i=0; i<taille; ++i)
            if (y[i]==k) nb++;
          if (nb==taille) ok++;
          if (MyGlobals::_Verbose>9)
            {
              res="0K"; if (nb!=taille) res="KO";
              cout<<res<<k<<" ";
            }
        }
      res="0K"; if (ok!=nb) res="MAUVAIS";
      if (MyGlobals::_Verbose>1) 
        cout<<"resultat "<<res<<" time(sec) "<<MPI_Wtime()-temps_debut<<endl;
      MPI_Request_free(&requete1);
    }
  //temps_fin=(MPI_WTIME()-temps_debut);
}

void MEDPARTITIONER::TestPersistantMpiRing(int taille, int nb)
{
  double temps_debut=MPI_Wtime();
  int befo, next, rang, wsize, tagbefo, tagnext;
  rang=MyGlobals::_Rank;
  wsize=MyGlobals::_World_Size;
  befo=rang-1; if (befo<0) befo=wsize-1;
  next=rang+1; if (next>=wsize) next=0;
  vector<int> x, y;
  tagbefo=111111+befo;
  tagnext=111111+rang;
  MPI_Request requete0, requete1;
  MPI_Status statut1, statut2;
  int ok=0;
  string res;
  //cout<<"ini|"<<rang<<'|'<<befo<<'|'<<next<<' ';
  {
    x.resize(taille);
    y.resize(taille);
    MPI_Ssend_init(&x[0], taille, MPI_INT, next, tagnext, MPI_COMM_WORLD , &requete0);
    MPI_Recv_init(&y[0], taille,  MPI_INT, befo, tagbefo, MPI_COMM_WORLD , &requete1);
    //cout<<"isr|"<<rang<<'|'<<requete0<<'|'<<requete1<<' ';
    for(int k=0; k<nb; k++)
      {
        for (int i=0; i<taille; ++i) x[i]=k+rang;
        //Envoi d’un gros message --> cela peut prendre du temps
        MPI_Start(&requete0);
        //Reception du gros message --> cela peut prendre du temps
        for (int i=0; i<taille; ++i) y[i]=(-1);
        MPI_Start(&requete1);
        //Traitement sequentiel independant de "x"
        //...
        //Traitement sequentiel independant de "y"
        //...
        //cout<<"dsr|"<<rang<<' ';
        MPI_Wait(&requete1, &statut1);
        //Traitement sequentiel dependant de "y"
        //...=f(y)
        int nb=0;
        for (int i=0; i<taille; ++i)
          if (y[i]==k+befo) nb++;
        if (nb==taille) ok++;
        if (MyGlobals::_Verbose>9)
          {
            res="0K"+IntToStr(rang); if (nb!=taille) res="KO"+IntToStr(rang);
            cout<<res<<k<<" ";
          }
        MPI_Wait(&requete0, &statut2);
        //Traitement sequentiel impliquant une modification de "x" en memoire
        //x=...
      }
    res="0K"; if (ok!=nb) res="MAUVAIS";
    temps_debut=MPI_Wtime()-temps_debut;
    MPI_Request_free(&requete1);
    MPI_Request_free(&requete0);
  }
  //temps_fin=(MPI_WTIME()-temps_debut);
  if (MyGlobals::_Verbose>1) 
    cout<<"resultat proc "<<rang<<" "<<res<<" time(sec) "<<temps_debut<<endl;
}
void MEDPARTITIONER::TestPersistantMpiRingOnCommSplit(int taille, int nb)
{
  double temps_debut=MPI_Wtime();
  int rang=MyGlobals::_Rank;
  MPI_Comm newcomm;
  int couleur=1;
  int rangMax=4;
  if (rang>=rangMax) couleur=MPI_UNDEFINED;
  cout<<"coul|"<<rang<<'|'<<couleur<<' ';
  //MPI_Comm_dup (MPI_COMM_WORLD, &newcomm) ;
  MPI_Comm_split(MPI_COMM_WORLD, couleur, rang, &newcomm);
  
  int befo, next, wsize, tagbefo, tagnext;
  wsize=rangMax;
  if (wsize>MyGlobals::_World_Size) wsize=MyGlobals::_World_Size;
  befo=rang-1; if (befo<0) befo=wsize-1;
  next=rang+1; if (next>=wsize) next=0;
  vector<int> x, y;
  tagbefo=111111+befo;
  tagnext=111111+rang;
  MPI_Request requete0, requete1;
  MPI_Status statut1, statut2;
  int ok=0;
  string res;
  
  //cout<<"ini|"<<rang<<'|'<<befo<<'|'<<next<<' ';
  if (couleur==1)
    {
      x.resize(taille);
      y.resize(taille);
      MPI_Ssend_init(&x[0], taille, MPI_INT, next, tagnext, newcomm , &requete0);
      MPI_Recv_init(&y[0], taille,  MPI_INT, befo, tagbefo, newcomm , &requete1);
      //cout<<"isr|"<<rang<<'|'<<requete0<<'|'<<requete1<<' ';
      for(int k=0; k<nb; k++)
        {
          for (int i=0; i<taille; ++i) x[i]=k+rang;
          //Envoi d’un gros message --> cela peut prendre du temps
          MPI_Start(&requete0);
          //Reception du gros message --> cela peut prendre du temps
          for (int i=0; i<taille; ++i) y[i]=(-1);
          MPI_Start(&requete1);
          //Traitement sequentiel independant de "x"
          //...
          //Traitement sequentiel independant de "y"
          //...
          //cout<<"dsr|"<<rang<<' ';
          MPI_Wait(&requete1, &statut1);
          //Traitement sequentiel dependant de "y"
          //...=f(y)
          int nb=0;
          for (int i=0; i<taille; ++i)
            if (y[i]==k+befo) nb++;
          if (nb==taille) ok++;
          if (MyGlobals::_Verbose>9)
            {
              res="0K"+IntToStr(rang); if (nb!=taille) res="KO"+IntToStr(rang);
              cout<<res<<k<<" ";
            }
          MPI_Wait(&requete0, &statut2);
          //Traitement sequentiel impliquant une modification de "x" en memoire
          //x=...
        }
      res="0K"; if (ok!=nb) res="MAUVAIS";
      temps_debut=MPI_Wtime()-temps_debut;
      MPI_Request_free(&requete1);
      MPI_Request_free(&requete0);
    }
  //MPI_Barrier(MPI_COMM_WORLD);
  cout<<"barrier|"<<rang<<"|"<<newcomm<<" ";
  if (couleur==1) MPI_Comm_free(&newcomm);
  //temps_fin=(MPI_WTIME()-temps_debut);
  if (MyGlobals::_Verbose>1) 
    cout<<"resultat proc "<<rang<<" "<<res<<" time(sec) "<<temps_debut<<endl;
}


