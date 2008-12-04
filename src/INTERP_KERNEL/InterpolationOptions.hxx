//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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
#ifndef _INTERPOLATIONOPTIONS_HXX_
#define _INTERPOLATIONOPTIONS_HXX_


namespace INTERP_KERNEL {
  typedef enum {Triangulation, Convex, Geometric2D, Generic} IntersectionType;
	/// Type describing the different ways in which the hexahedron can be split into tetrahedra.
	/// The PLANAR_* policies persume that each face is to be considered planar, while the general
	/// policies make no such hypothesis. The integer at the end gives the number of tetrahedra
	/// that result from the split.
	typedef enum  { PLANAR_FACE_5 = 5, PLANAR_FACE_6 = 6, GENERAL_24 = 24, GENERAL_48 = 48 } SplittingPolicy;

	
	class InterpolationOptions{
	private :
		int _printLevel ;
		IntersectionType _intersectionType;
		double _precision;
		double _medianPlane ;
		bool _doRotate ;
		double _boundingBoxAdjustment ;
		int _orientation ;
		SplittingPolicy _splittingPolicy ;

	public:
		InterpolationOptions()
		{init();
		}
		int getPrintLevel() const {return _printLevel;}
	  void setPrintLevel(int pl) {_printLevel=pl;}

		IntersectionType getIntersectionType() const{return InterpolationOptions::_intersectionType;}
		void setIntersectionType(IntersectionType it) {InterpolationOptions::_intersectionType=it;}

		double getPrecision() const {return InterpolationOptions::_precision;}
		void setPrecision(double p) {InterpolationOptions::_precision=p;}

		double getMedianPlane() {return InterpolationOptions::_medianPlane;}
		void setMedianPlane(double mp) {InterpolationOptions::_medianPlane=mp;}
		
		bool getDoRotate() {return InterpolationOptions::_doRotate;}
		void setDoRotate( bool dr) {InterpolationOptions::_doRotate = dr;}
		
		double getBoundingBoxAdjustment() {return InterpolationOptions::_boundingBoxAdjustment;}
		void setBoundingBoxAdjustment(double bba) {InterpolationOptions::_boundingBoxAdjustment=bba;}
		
		int getOrientation() {return InterpolationOptions::_orientation;}
		void setOrientation(int o) {InterpolationOptions::_orientation=o;}
		
		SplittingPolicy getSplittingPolicy() {return _splittingPolicy;}
		void setSplittingPolicy(SplittingPolicy sp) {_splittingPolicy=sp;}
		void init()
		{	
			 _printLevel=0;
			_intersectionType=Triangulation;
			_precision=1e-12;;
			_medianPlane =0.5;
			_doRotate =true;
			_boundingBoxAdjustment =0.1;
			_orientation =0;
			_splittingPolicy =GENERAL_48;
		}
	};

}
#endif
