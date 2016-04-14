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

#include "MEDPARTITIONER_Utils.hxx"

#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"
#include "CellModel.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "InterpKernelException.hxx"
#include "MCAuto.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace MEDPARTITIONER;

/*!
 * not optimized but suffisant
 * return empty vector if i am not target
 */
std::vector<std::string> MEDPARTITIONER::SendAndReceiveVectorOfString(const std::vector<std::string>& vec, const int source, const int target)
{
  int rank=MyGlobals::_Rank;

  MPI_Status status;
  int tag = 111001;
  if (rank == source)
    {
      std::string str=SerializeFromVectorOfString(vec);
      int size=str.length();
      MPI_Send( &size, 1, MPI_INT, target, tag, MPI_COMM_WORLD );
      MPI_Send( (void*)str.data(), str.length(), MPI_CHAR, target, tag+100, MPI_COMM_WORLD );
    }
  
  int recSize=0;
  if (rank == target)
    {
      MPI_Recv(&recSize, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
      std::string recData(recSize,'x');
      MPI_Recv((void*)recData.data(), recSize, MPI_CHAR, source, tag+100, MPI_COMM_WORLD, &status);
      return DeserializeToVectorOfString(recData); //not empty one for target proc
    }
  std::vector<std::string> res;
  return res; //empty one for other proc
}

/*!
 * strings NO need all same size!!!!
 */
std::vector<std::string> MEDPARTITIONER::AllgathervVectorOfString(const std::vector<std::string>& vec)
{
  if (MyGlobals::_World_Size==1) //nothing to do
    return vec;

  int world_size=MyGlobals::_World_Size;
  std::string str=SerializeFromVectorOfString(vec);
  
  std::vector<int> indexes(world_size);
  int size=str.length();
  MPI_Allgather(&size, 1, MPI_INT, 
                &indexes[0], 1, MPI_INT, MPI_COMM_WORLD);
  
  //calcul of displacement
  std::vector< int > disp(1,0);
  for (int i=0; i<world_size; i++) disp.push_back( disp.back() + indexes[i] );
  
  std::string recData(disp.back(),'x');
  MPI_Allgatherv((void*)str.data(), str.length(), MPI_CHAR,
                 (void*)recData.data(), &indexes[0], &disp[0], MPI_CHAR,
                 MPI_COMM_WORLD);
  
  //really extraordinary verbose for debug
  std::vector<std::string> deserial=DeserializeToVectorOfString(recData);
  if (MyGlobals::_Verbose>1000) 
    {
      std::cout << "proc "<<MyGlobals::_Rank<<" : receive '" << recData << "'" << std::endl;
      std::cout << "deserialize is : a vector of size " << deserial.size() << std::endl;
      std::cout << ReprVectorOfString(deserial) << std::endl;
    }
  return deserial;
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
    std::cout << "proc " << MyGlobals::_Rank << " : --> SendDoubleVec " << size << std::endl;
#ifdef HAVE_MPI
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
#ifdef HAVE_MPI
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000) 
    std::cout << "proc " << MyGlobals::_Rank << " : <-- RecvDoubleVec " << size << std::endl;
  std::vector<double>* vec=new std::vector<double>;
  vec->resize(size);
  MPI_Recv(&vec[0], size, MPI_DOUBLE, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return vec;
}

void MEDPARTITIONER::RecvDoubleVec(std::vector<double>& vec, const int source)
{
  int tag = 111002;
  int size;
#ifdef HAVE_MPI
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    std::cout<< "proc " << MyGlobals::_Rank << " : <-- RecvDoubleVec " << size << std::endl;;
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
    std::cout << "proc " << MyGlobals::_Rank << " : --> SendIntVec " << size << std::endl;
#ifdef HAVE_MPI
  MPI_Send(&size, 1, MPI_INT, target, tag, MPI_COMM_WORLD);
  MPI_Send(const_cast<int*>(&vec[0]), size,MPI_INT, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill vector<int> vec.
  To be used with \a SendIntVec method.
  \param vec vector that is filled
  \param source processor id of the incoming messages
*/
std::vector<int> *MEDPARTITIONER::RecvIntVec(const int source)
{
  int tag = 111003;
  int size;
#ifdef HAVE_MPI
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    std::cout << "proc " << MyGlobals::_Rank << " : <-- RecvIntVec " << size << std::endl;
  std::vector<int> *vec=new std::vector<int>;
  vec->resize(size);
  MPI_Recv(&vec[0], size, MPI_INT, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return vec;
}

void MEDPARTITIONER::RecvIntVec(std::vector<int>& vec, const int source)
{
  int tag = 111003;
  int size;
#ifdef HAVE_MPI
  MPI_Status status;  
  MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    std::cout << "proc " << MyGlobals::_Rank << " : <-- RecvIntVec " << size << std::endl;
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
void MEDPARTITIONER::SendDataArrayInt(const MEDCoupling::DataArrayInt *da, const int target)
{
  if (da==0)
    throw INTERP_KERNEL::Exception("Problem send DataArrayInt* NULL");
  int tag = 111004;
  int size[3];
  size[0]=da->getNbOfElems();
  size[1]=da->getNumberOfTuples();
  size[2]=da->getNumberOfComponents();
  if (MyGlobals::_Verbose>1000) 
    std::cout << "proc " << MyGlobals::_Rank << " : --> SendDataArrayInt " << size[0] << std::endl;
#ifdef HAVE_MPI
  MPI_Send(&size, 3, MPI_INT, target, tag, MPI_COMM_WORLD);
  const int *p=da->getConstPointer();
  MPI_Send(const_cast<int*>(&p[0]), size[0] ,MPI_INT, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill dataArrayInt da.
  To be used with \a SendIntVec method.
  \param da dataArrayInt that is filled
  \param source processor id of the incoming messages
*/
MEDCoupling::DataArrayInt *MEDPARTITIONER::RecvDataArrayInt(const int source)
{
  int tag = 111004;
  int size[3];
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Recv(size, 3, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    std::cout << "proc " << MyGlobals::_Rank << " : <-- RecvDataArrayInt " << size[0] << std::endl;
  if (size[0]!=(size[1]*size[2]))
    throw INTERP_KERNEL::Exception("Problem in RecvDataArrayInt incoherent sizes");
  MEDCoupling::DataArrayInt* da=MEDCoupling::DataArrayInt::New();
  da->alloc(size[1],size[2]);
  int *p=da->getPointer();
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
void MEDPARTITIONER::SendDataArrayDouble(const MEDCoupling::DataArrayDouble *da, const int target)
{
  if (da==0)
    throw INTERP_KERNEL::Exception("Problem send DataArrayDouble* NULL");
  int tag = 111005;
  int size[3];
  size[0]=da->getNbOfElems();
  size[1]=da->getNumberOfTuples();
  size[2]=da->getNumberOfComponents();
  if (MyGlobals::_Verbose>1000) 
    std::cout << "proc " << MyGlobals::_Rank << " : --> SendDataArrayDouble " << size[0] << std::endl;
#ifdef HAVE_MPI
  MPI_Send(&size, 3, MPI_INT, target, tag, MPI_COMM_WORLD);
  const double *p=da->getConstPointer();
  MPI_Send(const_cast<double*>(&p[0]), size[0] ,MPI_DOUBLE, target, tag+100, MPI_COMM_WORLD);
#endif
}

/*! Receives messages from proc \a source to fill dataArrayDouble da.
  To be used with \a SendDoubleVec method.
  \param da dataArrayDouble that is filled
  \param source processor id of the incoming messages
*/
MEDCoupling::DataArrayDouble* MEDPARTITIONER::RecvDataArrayDouble(const int source)
{
  int tag = 111005;
  int size[3];
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Recv(size, 3, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  if (MyGlobals::_Verbose>1000)
    std::cout << "proc " << MyGlobals::_Rank << " : <-- RecvDataArrayDouble " << size[0] << std::endl;
  if (size[0]!=(size[1]*size[2]))
    throw INTERP_KERNEL::Exception("Problem in RecvDataArrayDouble incoherent sizes");
  MEDCoupling::DataArrayDouble* da=MEDCoupling::DataArrayDouble::New();
  da->alloc(size[1],size[2]);
  double *p=da->getPointer();
  MPI_Recv(const_cast<double*>(&p[0]), size[0], MPI_DOUBLE, source, tag+100, MPI_COMM_WORLD, &status);
#endif
  return da;
}

void MEDPARTITIONER::TestVectorOfStringMpi()
{
  int rank=MyGlobals::_Rank;
  int world_size=MyGlobals::_World_Size;
  std::vector<std::string> myVector;
  std::ostringstream oss;
  oss << "hello from " << std::setw(5) << rank << " " << std::string(rank+1,'n') <<
    " next is an empty one";
  myVector.push_back(oss.str());
  myVector.push_back("");
  myVector.push_back("next is an singleton");
  myVector.push_back("1");
  
  if (rank==0)
    {
      std::string s0=SerializeFromVectorOfString(myVector);
      std::vector<std::string> res=DeserializeToVectorOfString(s0);
      if (res.size()!=myVector.size()) 
        throw INTERP_KERNEL::Exception("Problem in (de)serialise VectorOfString incoherent sizes");
      for (std::size_t i=0; i<myVector.size(); i++)
        if (res[i]!=myVector[i])
          throw INTERP_KERNEL::Exception("Problem in (de)serialise VectorOfString incoherent elements");
    }

  for (int i=0; i<world_size; i++)
    {
      for (int j=0; j<world_size; j++)
        {
          std::vector<std::string> res=SendAndReceiveVectorOfString(myVector, i, j);
          if ((rank==j) && MyGlobals::_Verbose>20)
            std::cout << "proc " << rank << " : receive \n" << ReprVectorOfString(res) << std::endl;
          if (rank==j)
            {
              if (res.size()!=myVector.size()) 
                throw INTERP_KERNEL::Exception("Problem in SendAndReceiveVectorOfString incoherent sizes");
              for (std::size_t ii=1; ii<myVector.size(); ii++) //first is different
                if (res[i]!=myVector[ii])
                  throw INTERP_KERNEL::Exception("Problem in SendAndReceiveVectorOfString incoherent elements");
            }
          else 
            {
              if (res.size()!=0) 
                throw INTERP_KERNEL::Exception("Problem in SendAndReceiveVectorOfString size have to be 0");
            }
        }
    }
  std::vector<std::string> res=AllgathervVectorOfString(myVector);
  //sometimes for test
  res=AllgathervVectorOfString(myVector);
  res=AllgathervVectorOfString(myVector);
  if (rank==0 && MyGlobals::_Verbose>20)
    std::cout << "proc " << rank << " : receive \n" << ReprVectorOfString(res) << std::endl;
  if (res.size()!=myVector.size()*world_size) 
    throw INTERP_KERNEL::Exception("Problem in AllgathervVectorOfString incoherent sizes");
  int jj=-1;
  for (int j=0; j<world_size; j++)
    {
      for (int i=0; i<(int)myVector.size(); i++)
        {
          jj=jj+1;
          if (i==0)
            continue; //first is different
          if (res[jj]!=myVector[i])
            throw INTERP_KERNEL::Exception("Problem in AllgathervVectorOfString incoherent elements");
        }
    }
  if (MyGlobals::_Verbose)
    std::cout << "proc " << rank << " : OK TestVectorOfStringMpi END" << std::endl;
}

void MEDPARTITIONER::TestMapOfStringIntMpi()
{
  int rank=MyGlobals::_Rank;
  std::map<std::string,int> myMap;
  myMap["one"]=1;
  myMap["two"]=22;  //a bug
  myMap["three"]=3;
  myMap["two"]=2; //last speaking override
  
  if (rank==0)
    {
      std::vector<std::string> v2=VectorizeFromMapOfStringInt(myMap);
      std::map<std::string,int> m3=DevectorizeToMapOfStringInt(v2);
      if (ReprMapOfStringInt(m3)!=ReprMapOfStringInt(myMap))
        throw INTERP_KERNEL::Exception("Problem in (de)vectorize MapOfStringInt");
    }
    
  std::vector<std::string> v2=AllgathervVectorOfString(VectorizeFromMapOfStringInt(myMap));
  if (rank==0 && MyGlobals::_Verbose>20)
    {
      std::cout << "v2 is : a vector of size " << v2.size() << std::endl;
      std::cout << ReprVectorOfString(v2) << std::endl;
      std::map<std::string,int> m2=DevectorizeToMapOfStringInt(v2);
      std::cout << "m2 is : a map of size " << m2.size() << std::endl;
      std::cout << ReprMapOfStringInt(m2) << std::endl;
    }
  if (MyGlobals::_Verbose)
    std::cout << "proc " << rank << " : OK TestMapOfStringIntMpi END" << std::endl;
}

void MEDPARTITIONER::TestMapOfStringVectorOfStringMpi()
{
  int rank=MyGlobals::_Rank;
  std::vector<std::string> myVector;
  std::ostringstream oss;
  oss << "hello from " << std::setw(5) << MyGlobals::_Rank << " " << std::string(rank+1,'n') << " next is an empty one";
  myVector.push_back(oss.str());
  myVector.push_back("");
  myVector.push_back("next is an singleton");
  myVector.push_back("1");
  
  if (rank==0)
    {
      std::map< std::string,std::vector<std::string> > m2;
      m2["first key"]=myVector;
      m2["second key"]=myVector;
      std::vector<std::string> v2=VectorizeFromMapOfStringVectorOfString(m2);
      std::map< std::string,std::vector<std::string> > m3=DevectorizeToMapOfStringVectorOfString(v2);
      if (rank==0 && MyGlobals::_Verbose>20)
        std::cout << "m2 is : a MapOfStringVectorOfString of size " << m2.size() << std::endl;
      std::cout << ReprMapOfStringVectorOfString(m2) << std::endl;
      std::cout << "v2 is : a vector of size " << v2.size() << std::endl;
      std::cout << ReprVectorOfString(v2) << std::endl;
      std::cout << "m3 is : a map of size "<<m3.size() << std::endl;
      std::cout << ReprMapOfStringVectorOfString(m3) << std::endl;
      if (ReprMapOfStringVectorOfString(m3)!=ReprMapOfStringVectorOfString(m2))
        throw INTERP_KERNEL::Exception("Problem in (de)vectorize MapOfStringVectorOfString");
    }
    
  std::map< std::string,std::vector<std::string> > m4;
  m4["1rst key"]=myVector;
  m4["2snd key"]=myVector;
  std::vector<std::string> v4=AllgathervVectorOfString(VectorizeFromMapOfStringVectorOfString(m4));
  if (rank==0 && MyGlobals::_Verbose>20)
    {
      std::map< std::string,std::vector<std::string> > m5=DevectorizeToMapOfStringVectorOfString(v4);
      std::map< std::string,std::vector<std::string> > m6=DeleteDuplicatesInMapOfStringVectorOfString(m5);
      std::cout<< "m5 is : a map of size "<<m5.size() << std::endl;
      std::cout<< ReprMapOfStringVectorOfString(m5) << std::endl;
      std::cout<< "m6 is : a map from m5 with deleteDuplicates of size " << m6.size() << std::endl;
      std::cout<< ReprMapOfStringVectorOfString(m6) << std::endl;
    }
  if (MyGlobals::_Verbose)
    std::cout<<"proc " << rank << " : OK TestMapOfStringVectorOfStringMpi END" << std::endl;
}

void MEDPARTITIONER::TestDataArrayMpi()
{
  int rank=MyGlobals::_Rank;
  //int
  {
    MEDCoupling::DataArrayInt* send=MEDCoupling::DataArrayInt::New();
    MEDCoupling::DataArrayInt* recv=0;
    int nbOfTuples=5;
    int numberOfComponents=3;
    send->alloc(nbOfTuples,numberOfComponents);
    std::vector<int> vals;
    for (int j=0; j<nbOfTuples; j++)
      for (int i=0; i<numberOfComponents; i++) vals.push_back((j+1)*10+i+1);
    std::copy(vals.begin(),vals.end(),send->getPointer());
    if (rank==0)
      SendDataArrayInt(send, 1);
    if (rank==1)
      recv=RecvDataArrayInt(0);
    if (rank==1 && MyGlobals::_Verbose>20)
      {
        std::cout << send->repr() << std::endl;
        std::cout << recv->repr() << std::endl;
      }
    if (rank==1)
      {
        if (send->repr()!=recv->repr())
          throw INTERP_KERNEL::Exception("Problem in send&recv DataArrayInt");
      }
    send->decrRef();
    if (rank==1)
      recv->decrRef();
  }
  //double
  {
    MEDCoupling::DataArrayDouble* send=MEDCoupling::DataArrayDouble::New();
    MEDCoupling::DataArrayDouble* recv=0;
    int nbOfTuples=5;
    int numberOfComponents=3;
    send->alloc(nbOfTuples,numberOfComponents);
    std::vector<double> vals;
    for (int j=0; j<nbOfTuples; j++)
      for (int i=0; i<numberOfComponents; i++) vals.push_back(double(j+1)+double(i+1)/10);
    std::copy(vals.begin(),vals.end(),send->getPointer());
    if (rank==0) SendDataArrayDouble(send, 1);
    if (rank==1) recv=RecvDataArrayDouble(0);
    if (rank==1 && MyGlobals::_Verbose>20)
      {
        std::cout << send->repr() << std::endl;
        std::cout << recv->repr() << std::endl;
      }
    if (rank==1)
      {
        if (send->repr()!=recv->repr())
          throw INTERP_KERNEL::Exception("Problem in send&recv DataArrayDouble");
      }
    send->decrRef();
    if (rank==1) recv->decrRef();
  }
  
  if (MyGlobals::_Verbose)
    std::cout << "proc " << rank << " : OK TestDataArrayMpi END" << std::endl;
}

void MEDPARTITIONER::TestPersistantMpi0To1(int taille, int nb)
{
  double temps_debut=MPI_Wtime();
  int rank=MyGlobals::_Rank;
  std::vector<int> x, y;
  int tag=111111;
  MPI_Request requete0, requete1;
  MPI_Status statut;
  int ok=0;
  std::string res;
  if (rank==0)
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
  else if (rank == 1)
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
          int nbb=0;
          for (int i=0; i<taille; ++i)
            if (y[i]==k)
              nbb++;
          if (nbb==taille)
            ok++;
          if (MyGlobals::_Verbose>9)
            {
              res="0K";
              if (nbb!=taille)
                res="KO";
              std::cout << res << k << " ";
            }
        }
      res="0K";
      if (ok!=nb)
        res="BAD";
      if (MyGlobals::_Verbose>1) 
        std::cout << "result " << res << " time(sec) " << MPI_Wtime()-temps_debut << std::endl;
      MPI_Request_free(&requete1);
    }
  //end_time=(MPI_WTIME()-start_time);
}

void MEDPARTITIONER::TestPersistantMpiRing(int taille, int nb)
{
  double temps_debut=MPI_Wtime();
  int befo, next, rank, wsize, tagbefo, tagnext;
  rank=MyGlobals::_Rank;
  wsize=MyGlobals::_World_Size;
  befo=rank-1; if (befo<0) befo=wsize-1;
  next=rank+1; if (next>=wsize) next=0;
  std::vector<int> x, y;
  tagbefo=111111+befo;
  tagnext=111111+rank;
  MPI_Request requete0, requete1;
  MPI_Status statut1, statut2;
  int ok=0;
  std::string res;
  //cout<<"ini|"<<rank<<'|'<<befo<<'|'<<next<<' ';
  {
    x.resize(taille);
    y.resize(taille);
    MPI_Ssend_init(&x[0], taille, MPI_INT, next, tagnext, MPI_COMM_WORLD , &requete0);
    MPI_Recv_init(&y[0], taille,  MPI_INT, befo, tagbefo, MPI_COMM_WORLD , &requete1);
    for(int k=0; k<nb; k++)
      {
        for (int i=0; i<taille; ++i) x[i]=k+rank;
        //Envoi d’un gros message --> cela peut prendre du temps
        MPI_Start(&requete0);
        //Reception du gros message --> cela peut prendre du temps
        for (int i=0; i<taille; ++i) y[i]=(-1);
        MPI_Start(&requete1);
        //Traitement sequentiel independant de "x"
        //...
        //Traitement sequentiel independant de "y"
        //...
        MPI_Wait(&requete1, &statut1);
        //Traitement sequentiel dependant de "y"
        //...=f(y)
        int nbb=0;
        for (int i=0; i<taille; ++i)
          if (y[i]==k+befo)
            nbb++;
        if (nbb==taille)
          ok++;
        if (MyGlobals::_Verbose>9)
          {
            res="0K"+IntToStr(rank);
            if (nbb!=taille)
              res="KO"+IntToStr(rank);
            std::cout << res << k << " ";
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
  //end_time=(MPI_WTIME()-start_time);
  if (MyGlobals::_Verbose>1) 
    std::cout << "result on proc " << rank << " " << res << " time(sec) " << temps_debut << std::endl;
}

void MEDPARTITIONER::TestPersistantMpiRingOnCommSplit(int size, int nb)
{
  double temps_debut=MPI_Wtime();
  int rank=MyGlobals::_Rank;
  MPI_Comm newcomm;
  int color=1;
  int rankMax=4;
  if (rank>=rankMax)
    color=MPI_UNDEFINED;
  //MPI_Comm_dup (MPI_COMM_WORLD, &newcomm) ;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &newcomm);
  
  int befo, next, wsize, tagbefo, tagnext;
  wsize=rankMax;
  if (wsize>MyGlobals::_World_Size)
    wsize=MyGlobals::_World_Size;
  befo=rank-1;
  if (befo<0)
    befo=wsize-1;
  next=rank+1;
  if (next>=wsize)
    next=0;
  std::vector<int> x, y;
  tagbefo=111111+befo;
  tagnext=111111+rank;
  MPI_Request requete0, requete1;
  MPI_Status statut1, statut2;
  int ok=0;
  std::string res;
  
  if (color==1)
    {
      x.resize(size);
      y.resize(size);
      MPI_Ssend_init(&x[0], size, MPI_INT, next, tagnext, newcomm , &requete0);
      MPI_Recv_init(&y[0], size,  MPI_INT, befo, tagbefo, newcomm , &requete1);
      for(int k=0; k<nb; k++)
        {
          for (int i=0; i<size; ++i)
            x[i]=k+rank;
          //Send of big message --> time consuming
          MPI_Start(&requete0);
          //Reception of big message --> time consuming
          for (int i=0; i<size; ++i)
            y[i]=-1;
          MPI_Start(&requete1);
          //Traitement sequentiel independant de "x"
          //...
          //Traitement sequentiel independant de "y"
          //...
          //cout<<"dsr|"<<rank<<' ';
          MPI_Wait(&requete1, &statut1);
          //Traitement sequentiel dependant de "y"
          //...=f(y)
          int nbb=0;
          for (int i=0; i<size; ++i)
            if (y[i]==k+befo)
              nbb++;
          if (nbb==size)
            ok++;
          if (MyGlobals::_Verbose>9)
            {
              res="0K"+IntToStr(rank);
              if (nbb!=size)
                res="KO"+IntToStr(rank);
              std::cout << res << k << " ";
            }
          MPI_Wait(&requete0, &statut2);
          //Traitement sequentiel impliquant une modification de "x" en memoire
          //x=...
        }
      res="0K";
      if (ok!=nb)
        res="MAUVAIS";
      temps_debut=MPI_Wtime()-temps_debut;
      MPI_Request_free(&requete1);
      MPI_Request_free(&requete0);
    }
  //MPI_Barrier(MPI_COMM_WORLD);
  if (color==1)
    MPI_Comm_free(&newcomm);
  if (MyGlobals::_Verbose>1) 
    std::cout << "resultat proc " << rank <<" " << res << " time(sec) " << temps_debut << std::endl;
}
