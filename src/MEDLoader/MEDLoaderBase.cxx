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

#include "MEDLoaderBase.hxx"
#include "InterpKernelException.hxx"

#include <sstream>
#include <fstream>
#include <cstring>
#include <iostream>

const char MEDLoaderBase::WHITE_SPACES[]=" \n";

int MEDLoaderBase::getStatusOfFile(const char *fileName)
{
  std::ifstream ifs;
  ifs.open(fileName);
  unsigned int res=0;
  if((ifs.rdstate() & std::ifstream::failbit)!=0)
    {
      res+=1;
      ifs.close();
    }
  std::ofstream ofs(fileName,std::ios_base::app);
  if((ofs.rdstate() & std::ofstream::failbit)!=0)
    {
      ofs.close();
      res+=2;
    }
  switch(res)
    {
    case 0:
      return EXIST_RW;
    case 1:
      {
        std::ifstream ifs2;
        ifs2.open(fileName);
        if((ifs2.rdstate() & std::ifstream::failbit)!=0)
          return EXIST_WRONLY;
        else
          return NOT_EXIST;
      }
    case 2:
      return EXIST_RDONLY;
    case 3:
      return DIR_LOCKED;
    default:
      throw INTERP_KERNEL::Exception("Internal error !");
    }
}

char *MEDLoaderBase::buildEmptyString(int lgth)
{
  char *ret=new char[lgth+1];
  std::fill(ret,ret+lgth,' ');
  ret[lgth]='\0';
  return ret;
}

std::string MEDLoaderBase::buildUnionUnit(const char *name, const char *unit)
{
  std::string ret(name);
  if(unit[0]=='\0')
    return ret;
  ret+=" (";
  ret+=unit;
  ret+=")";
  return ret;
}

void MEDLoaderBase::splitIntoNameAndUnit(const std::string& s, std::string& name, std::string& unit)
{
  std::string::size_type f1=s.find_first_of('(');
  std::string::size_type f2=s.find_last_of(')');
  if(f1!=std::string::npos && f2!=std::string::npos)
    {
      if(f1<f2)
        {
          name=s.substr(0,f1);
          unit=s.substr(f1+1,f2-f1-1);
          strip(name);
          strip(unit);
          return;
        }
    }
  name=s;
  unit="";
}

void MEDLoaderBase::strip(std::string& s)
{
  std::string::size_type f1=s.find_first_not_of(' ');
  std::string::size_type f2=s.find_last_not_of(' ');
  s=s.substr(f1,f2-f1+1);
}

/*!
 * This method operates a safe copy from 'src' to 'dest' by checking the size of 'src' before trying to copy.
 * If size of 'src' string is higher than 'maxLgth' the behaviour is dependant from 'behaviour' parameter.
 * If 'behaviour' equals 0 an exception is thrown. If 'behaviour' equals 1 an attempt of zipping of string will be done
 * ( see zipString to have more details).
 */
void MEDLoaderBase::safeStrCpy(const char *src, int maxLgth, char *dest, int behaviour) throw(INTERP_KERNEL::Exception)
{
  if(strlen(src)>maxLgth)
    {
      if(behaviour==0 || behaviour>1)
        {
          std::ostringstream oss; oss << "A string : \"" << src << "\" has been detected to be too long for MED File ( > " << maxLgth << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      else if(behaviour==1)
        {
          std::string s=zipString(src,maxLgth);
          std::cerr << "A string : \"" << src << "\" has been detected to be too long for MED File ( > " << maxLgth << ") : ";
          std::cerr << "zipping to : " << s << "\n";
          strcpy(dest,s.c_str());
          return ;
        }
    }
  strcpy(dest,src);
}

std::string MEDLoaderBase::buildStringFromFortran(const char *expr, int lgth)
{
  std::string ret(expr,lgth);
  std::string whiteSpaces(WHITE_SPACES);
  std::size_t lgthReal=strlen(ret.c_str());
  std::string ret2=ret.substr(0,lgthReal);
  std::size_t found=ret2.find_last_not_of(whiteSpaces);
  if (found!=std::string::npos)
    ret2.erase(found+1);
  else
    ret2.clear();//ret is all whitespace
  return ret2;
}

/*!
 * This method given the target size to respect 'sizeToRespect' tries to reduce size of 'src' string.
 * This method uses several soft methods to do its job. But if it fails a simple cut of string will be performed.
 */
std::string MEDLoaderBase::zipString(const char *src, int sizeToRespect)
{
  std::string s(src);
  strip(s);
  if(s.length()<=sizeToRespect)
    return s;
  s=src;
  zipEqualConsChar(s,3);
  if(s.length()<=sizeToRespect)
    return s;
  s=src;
  zipEqualConsChar(s,2);
  if(s.length()<=sizeToRespect)
    return s;
  s=src;
  return s.substr(0,sizeToRespect);
}

/*!
 * This method see if there is in 's' more than 'minConsSmChar' consecutive same character.
 * If yes, the group will be zipped using only one character for this group.
 * If no such group is found, s remains unchanged.
 */
void MEDLoaderBase::zipEqualConsChar(std::string& s, int minConsSmChar)
{
  for(std::string::iterator it=s.begin();it!=s.end();it++)
    {
      char tmp=*it;
      int sz=1;
      for(std::string::iterator it2=it+1;it2!=s.end() && *it2==tmp;it2++)
        sz++;
      if(sz>=minConsSmChar)
        s.erase(it+1,it+sz);
    }
}

