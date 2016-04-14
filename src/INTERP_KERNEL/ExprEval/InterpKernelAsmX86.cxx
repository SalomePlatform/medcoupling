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
// Author : Anthony Geay (CEA/DEN)

#include "InterpKernelAsmX86.hxx"

#include <cstring>
#include <sstream>
#include <algorithm>

#ifdef _POSIX_MAPPED_FILES
#include <sys/mman.h>
#else
#ifdef WIN32
#include <windows.h>
#endif
#endif

const char *INTERP_KERNEL::AsmX86::OPS[NB_OF_OPS]={"mov","push","pop","fld","faddp","fsubp","fmulp","fdivp","fcos","fsin","fabs","fchs","fsqrt","sub","add","ret","leave","movsd","fst"};

std::vector<char> INTERP_KERNEL::AsmX86::convertIntoMachineLangage(const std::vector<std::string>& asmb) const
{
  std::vector<char> ret;
  for(std::vector<std::string>::const_iterator iter=asmb.begin();iter!=asmb.end();iter++)
    convertOneInstructionInML(*iter,ret);
  return ret;
}

char *INTERP_KERNEL::AsmX86::copyToExecMemZone(const std::vector<char>& ml, unsigned& offset) const
{
  char *ret=0;
  int lgth=ml.size();
#ifdef _POSIX_MAPPED_FILES
# ifdef __APPLE__
  ret=(char *)mmap(0,lgth,PROT_EXEC | PROT_WRITE,MAP_ANON | MAP_PRIVATE,-1,0);
# else
  ret=(char *)mmap(0,lgth,PROT_EXEC | PROT_WRITE,MAP_ANONYMOUS | MAP_PRIVATE,-1,0);
# endif
#else
#ifdef WIN32
  HANDLE h=CreateFileMapping(INVALID_HANDLE_VALUE,NULL,PAGE_EXECUTE_READWRITE,0,lgth,NULL);
  ret=(char *)MapViewOfFile(h,FILE_MAP_EXECUTE | FILE_MAP_READ | FILE_MAP_WRITE,0,0,lgth);
#endif
#endif
  if(ret)
    std::copy(ml.begin(),ml.end(),ret);
  return ret;
}

void INTERP_KERNEL::AsmX86::convertOneInstructionInML(const std::string& inst, std::vector<char>& ml) const
{
  std::string::size_type pos=inst.find_first_of(' ');
  std::string op;
  std::string param;
  if(pos!=std::string::npos)
    {
      op=inst.substr(0,pos);
      param=inst.substr(pos+1);
    }
  else
    op=inst;
  int id=0;
  for(const char **it=OPS;it!=OPS+NB_OF_OPS;it++,id++)
    {
      std::string tmp(*it);
      if(op==tmp)
        break;
    }
  switch(id)
    {
    case 0:
      convertMov(param,ml);
      break;
    case 1:
      convertPush(param,ml);
      break;
    case 2:
      convertPop(param,ml);
      break;
    case 3:
      convertFld(param,ml);
      break;
    case 4:
      convertFaddp(param,ml);
      break;
    case 5:
      convertFsubp(param,ml);
      break;
    case 6:
      convertFmulp(param,ml);
      break;
    case 7:
      convertFdivp(param,ml);
      break;
    case 8:
      convertFcos(param,ml);
      break;
    case 9:
      convertFsin(param,ml);
      break;
    case 10:
      convertFabs(param,ml);
      break;
    case 11:
      convertFchs(param,ml);
      break;
    case 12:
      convertFsqrt(param,ml);
      break;
    case 13:
      convertSub(param,ml);
      break;
    case 14:
      convertAdd(param,ml);
      break;
    case 15:
      convertRet(param,ml);
      break;
    case 16:
      convertLeave(param,ml);
      break;
    case 17:
      convertMovsd(param,ml);
      break;
    case 18:
      convertFst(param,ml);
      break;
    default:
      {
        std::ostringstream oss; oss << "Unrecognized op : " << op << " in assembly line : " << inst;
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
}

#include <iostream>

void INTERP_KERNEL::AsmX86::convertMov(const std::string& inst, std::vector<char>& ml)
{
  const char ASM1[]="ebp,esp";
  const unsigned char ML1[2]={0x89,0xe5};
  if(inst==ASM1)
    {
      ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
      return ;
    }
  const char ASM2[]="rbp,rsp";
  const unsigned char ML2[3]={0x48,0x89,0xe5};
  if(inst==ASM2)
    {
      ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
      return ;
    }
  std::string::size_type pos=inst.find_first_of(' ');
  if(pos==std::string::npos)
    {
      std::ostringstream oss; oss << "not recognized instruction mov : " << inst;
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::string inst2=inst.substr(pos+1);
  pos=inst2.find_first_of(',');
  if(pos==std::string::npos)
    {
      std::ostringstream oss; oss << "not recognized instruction mov : " << inst;
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::string inst3=inst2.substr(0,pos);
  std::string inst4=inst2.substr(pos+1);
  convertMovToEsp(inst3,inst4,ml);
}

void INTERP_KERNEL::AsmX86::convertMovToEsp(const std::string& inst1, const std::string& inst2, std::vector<char>& ml)
{
  if(inst1[0]!='[' || inst1[inst1.length()-1]!=']')
    throw INTERP_KERNEL::Exception("not recognized convertMovToEsp exp !");
  std::string inst1bis=inst1.substr(1,inst1.length()-2);
  const char ASM1[]="esp";
  const unsigned char ML1[3]={0xc7,0x04,0x24};
  if(inst1bis==ASM1)
    {//mov dword [esp],0x3ff3c0ca
      ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
      appendAddress(inst2,4,ml);
      return ;
    }
  if(inst1bis.substr(0,3)==ASM1)
    {
      if(inst1bis[3]=='+')
        {//mov dword [esp+4],0x3ff3c0ca
          const unsigned char ML2[3]={0xc7,0x44,0x24};
          ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
          std::string::size_type pos=inst1bis.find_first_of(']');
          std::string inst1_1=inst1bis.substr(4,pos-4-1);
          appendAddress(inst1_1,1,ml);
          appendAddress(inst2,4,ml);
          return;
        }
      else
        throw INTERP_KERNEL::Exception("Not recognized exp : mov [esp@..],...");
    }
  const char ASM3[]="rsp";
  const unsigned char ML3[3]={0xc7,0x04,0x24};
  if(inst1bis==ASM3)
    {//mov dword [rsp],0x3ff3c0ca
      ml.insert(ml.end(),ML3,ML3+sizeof(ML3));
      appendAddress(inst2,4,ml);
      return ;
    }
  if(inst1bis.substr(0,3)==ASM3)
    {
      if(inst1bis[3]=='+')
        {//mov dword [rsp+4],0x3ff3c0ca
          const unsigned char ML2[3]={0xc7,0x44,0x24};
          ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
          std::string::size_type pos=inst1bis.find_first_of(']');
          std::string inst1_1=inst1bis.substr(4,pos-4-1);
          appendAddress(inst1_1,1,ml);
          appendAddress(inst2,4,ml);
          return;
        }
      else
        throw INTERP_KERNEL::Exception("Not recognized exp : mov [esp@..],...");
    }
  throw INTERP_KERNEL::Exception("Not recognized exp : mov");
}

void INTERP_KERNEL::AsmX86::convertPush(const std::string& inst, std::vector<char>& ml)
{
  std::string::size_type pos=inst.find_first_of(' ');
  std::string inst2=inst.substr(pos+1);
  const char ASM1[]="ebp";
  const unsigned char ML1[1]={0x55};
  if(inst2==ASM1)
    {//push ebp
      ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
      return ;
    }
  const char ASM2[]="ebx";
  const unsigned char ML2[1]={0x53};
  if(inst2==ASM2)
    {//push ebx
      ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
      return ;
    }
  const char ASM3[]="rbp";
  const unsigned char ML3[1]={0x55};
  if(inst2==ASM3)
    {//push rbp
      ml.insert(ml.end(),ML3,ML3+sizeof(ML3));
      return ;
    }
  throw INTERP_KERNEL::Exception("Unrecognized push instruction");
}

void INTERP_KERNEL::AsmX86::convertPop(const std::string& inst, std::vector<char>& ml)
{
  std::string::size_type pos=inst.find_first_of(' ');
  std::string inst2=inst.substr(pos+1);
  const char ASM1[]="ebp";
  const unsigned char ML1[1]={0x5d};
  if(inst2==ASM1)
    {//push ebp
      ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
      return ;
    }
  const char ASM2[]="ebx";
  const unsigned char ML2[1]={0x5b};
  if(inst2==ASM2)
    {//push ebx
      ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
      return ;
    }
  throw INTERP_KERNEL::Exception("Unrecognized pop instruction");
}

void INTERP_KERNEL::AsmX86::convertFld(const std::string& inst, std::vector<char>& ml)
{
  std::string::size_type pos=inst.find_first_of(' ');
  std::string params=inst.substr(pos+1);
  std::string params2=params.substr(1,params.length()-2);
  if(params2.substr(0,3)=="esp")
    {
      const unsigned char ML1[3]={0xdd,0x04,0x24};
      if(params2.length()==3)
        {//fld qword [esp]
          ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
          return ;
        }
      pos=params2.find_first_of('+');
      if(pos!=std::string::npos)
        {//fld qword [esp+@]
          ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
          std::string params3=params2.substr(pos+1);
          appendAddress(params3,1,ml);
          return ;
        }
      throw INTERP_KERNEL::Exception("Unrecognized fld esp...");
    }
  if(params2.substr(0,3)=="ebp")
    {
      const unsigned char ML2[2]={0xdd,0x45};
      if(params2.length()==3)
        {//fld qword [ebp]
          ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
          ml.push_back(0);
          return ;
        }
      pos=params2.find_first_of('+');
      if(pos!=std::string::npos)
        {//fld qword [esp+@]
          ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
          std::string params3=params2.substr(pos+1);
          appendAddress(params3,1,ml);
          return ;
        }
      throw INTERP_KERNEL::Exception("Unrecognized fld ebp...");
    }
  if(params2.substr(0,3)=="rsp")
    {
      const unsigned char ML2[3]={0xdd,0x04,0x24};
      ml.insert(ml.end(),ML2,ML2+sizeof(ML2));// to improve ! no fully managed !
      return ;
    }
  throw INTERP_KERNEL::Exception("Unrecognized fld instruction");
}

void INTERP_KERNEL::AsmX86::convertFaddp(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML1[2]={0xde,0xc1};
  ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
}

void INTERP_KERNEL::AsmX86::convertFsubp(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML1[2]={0xde,0xe9};
  ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
}

void INTERP_KERNEL::AsmX86::convertFmulp(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML1[2]={0xde,0xc9};
  ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
}

void INTERP_KERNEL::AsmX86::convertFdivp(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML1[2]={0xde,0xf9};
  ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
}

void INTERP_KERNEL::AsmX86::convertFcos(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[2]={0xd9,0xff};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertFsin(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[2]={0xd9,0xfe};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertFabs(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[2]={0xd9,0xe1};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertFchs(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[2]={0xd9,0xe0};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertFsqrt(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[2]={0xd9,0xfa};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertSub(const std::string& inst, std::vector<char>& ml)
{
  if(inst.substr(0,4)=="esp,")
    {
      const unsigned char ML[2]={0x81,0xec};
      ml.insert(ml.end(),ML,ML+sizeof(ML));
      std::string inst2=inst.substr(4);
      appendAddress(inst2,4,ml);
      return;
    }
  if(inst.substr(0,4)=="rsp,")
    {
      const unsigned char ML[4]={0x48,0x83,0xec,0x08};
      ml.insert(ml.end(),ML,ML+sizeof(ML)); // to improve 8 statically put (last of element of ML) !!!!
      return;
    }
  throw INTERP_KERNEL::Exception("Not recognized sub instruction.");
}

void INTERP_KERNEL::AsmX86::convertAdd(const std::string& inst, std::vector<char>& ml)
{
  if(inst.substr(0,4)=="esp,")
    {
      const unsigned char ML[2]={0x81,0xc4};
      ml.insert(ml.end(),ML,ML+sizeof(ML));
      std::string inst2=inst.substr(4);
      appendAddress(inst2,4,ml);
      return;
    }
  if(inst.substr(0,4)=="rsp,")
    {
      const unsigned char ML[4]={0x48,0x83,0xc4,0x08};
      ml.insert(ml.end(),ML,ML+sizeof(ML)); // to improve 8 statically put (last of element of ML) !!!!
      return;
    }
  throw INTERP_KERNEL::Exception("Not recognized add instruction.");
}

void INTERP_KERNEL::AsmX86::convertRet(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[1]={0xc3};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertLeave(const std::string& inst, std::vector<char>& ml)
{
  const unsigned char ML[1]={0xc9};
  ml.insert(ml.end(),ML,ML+sizeof(ML));
}

void INTERP_KERNEL::AsmX86::convertMovsd(const std::string& inst, std::vector<char>& ml)
{
  const char ASM1[]="[rsp],xmm0";
  const unsigned char ML1[5]={0xf2,0x0f,0x11,0x04,0x24};
  if(inst==ASM1)
    {
      ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
      return ;
    }
  const char ASM2[]="xmm0,[rsp]";
  const unsigned char ML2[5]={0xf2,0x0f,0x10,0x04,0x24};
  if(inst==ASM2)
    {
      ml.insert(ml.end(),ML2,ML2+sizeof(ML2));
      return ;
    }
  std::ostringstream oss; oss << "not recognized instruction movsd : " << inst;
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

void INTERP_KERNEL::AsmX86::convertFst(const std::string& inst, std::vector<char>& ml)
{
  const char ASM1[]="qword [rsp]";
  const unsigned char ML1[3]={0xdd,0x14,0x24};
  if(inst==ASM1)
    {
      ml.insert(ml.end(),ML1,ML1+sizeof(ML1));
      return ;
    }
  std::ostringstream oss; oss << "not recognized instruction fst : " << inst;
  throw INTERP_KERNEL::Exception(oss.str().c_str());
  //tony
}


void INTERP_KERNEL::AsmX86::appendAddress(const std::string& addr, int nbOfByte, std::vector<char>& ml)
{
  int i,j;
  char v;
  std::istringstream iss(addr);
  if(addr.length()>2)
    {
      if(addr[0]=='0' && addr[1]=='x')
        iss >> std::hex;
    }
  iss >> i;
  for(int k=0;k<nbOfByte;k++)
    {
      j=i&255;
      v=j;
      ml.push_back(v);
      i>>=8;
    }
}
