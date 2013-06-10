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

#include "MEDMeshMaker.hxx"

#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Meshing.hxx"

MEDMEM::MESH* MEDMeshMaker(int dim, int nbedge, MED_EN::medGeometryElement type)
{
  MEDMEM::MESHING* mesh=new MEDMEM::MESHING();
  int nbnodes;
  int nbelems;
  switch (dim)
    {
    case 2: 
      nbnodes=(nbedge+1)*(nbedge+1);
      if(type==MED_EN::MED_QUAD4)
        nbelems=(nbedge*nbedge);
      else
        throw MEDMEM::MEDEXCEPTION("MEDMeshMaker: type not impletmented");
      break;
    case 3:
      nbnodes=(nbedge+1)*(nbedge+1)*(nbedge+1);
      if (type==MED_EN::MED_HEXA8)
        nbelems= nbedge*nbedge*nbedge;
      else
        throw MEDMEM::MEDEXCEPTION("MEDMeshMaker: type not impletmented");
      break;
    }
  double* coords = new double[dim*nbnodes];
  int nz;
  if (dim==2) nz =1; else nz=nbedge+1;
  {
    for (int ix=0; ix < nbedge+1; ix++)
      for (int iy=0; iy<nbedge+1; iy++)
        for (int iz=0; iz<nz;iz++)
          {
            int inode=(ix*(nbedge+1)*nz+iy*nz+iz);
            coords[inode*dim]=double(ix)/double(nbedge);
            coords[inode*dim+1]=double(iy)/double(nbedge);
            if (dim==3)
              coords[inode*dim+2]=double(iz)/double(nbedge);
          }
  }
  mesh->setCoordinates(dim, nbnodes,coords,"CARTESIAN",MED_EN::MED_FULL_INTERLACE);
  delete [] coords;
  mesh->setNumberOfTypes(1,MED_EN::MED_CELL);
  mesh->setTypes(&type,MED_EN::MED_CELL);
  mesh->setNumberOfElements(&nbelems,MED_EN::MED_CELL);
  
  int* conn = new int [nbelems*(type%100)];
  if (dim==2)
    {
      for (int ix=0; ix<nbedge; ix++)
        for (int iy=0; iy<nbedge; iy++)
          {
            int ielem=(ix*nbedge+iy);
            conn [ielem*4]=ix*(nbedge+1)+iy+1;
            conn [ielem*4+1]=ix*(nbedge+1)+iy+1+1;
            conn [ielem*4+2]=(ix+1)*(nbedge+1)+iy+1+1;
            conn [ielem*4+3]=(ix+1)*(nbedge+1)+iy+1;                               
          }
    }
  if (dim==3)
    {
      for (int ix=0; ix<nbedge; ix++)
        for (int iy=0; iy<nbedge; iy++)
          for (int iz=0; iz<nbedge; iz++)
            {
              int ielem=(ix*nbedge*nbedge+iy*nbedge+iz);
              conn [ielem*8]=ix*(nbedge+1)*(nbedge+1)+iy*(nbedge+1)+iz+1;
              conn [ielem*8+1]=(ix+1)*(nbedge+1)*(nbedge+1)+iy*(nbedge+1)+iz+1;
              conn [ielem*8+2]=(ix+1)*(nbedge+1)*(nbedge+1)+(iy+1)*(nbedge+1)+iz+1;
              conn [ielem*8+3]=ix*(nbedge+1)*(nbedge+1)+(iy+1)*(nbedge+1)+iz+1;
              conn [ielem*8+4]=ix*(nbedge+1)*(nbedge+1)+iy*(nbedge+1)+iz+1+1;
              conn [ielem*8+5]=(ix+1)*(nbedge+1)*(nbedge+1)+iy*(nbedge+1)+iz+1+1;
              conn [ielem*8+6]=(ix+1)*(nbedge+1)*(nbedge+1)+(iy+1)*(nbedge+1)+iz+1+1;
              conn [ielem*8+7]=ix*(nbedge+1)*(nbedge+1)+(iy+1)*(nbedge+1)+iz+1+1;
            }
    }
  mesh->setConnectivity(MED_EN::MED_CELL,type,conn);
  delete [] conn;
  return mesh;
}
