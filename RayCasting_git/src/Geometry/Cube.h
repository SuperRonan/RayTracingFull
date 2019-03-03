#ifndef _Geometry_Cube_H
#define _Geometry_Cube_H

#include <Geometry/Geometry.h>
#include <Geometry/Square.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Cube
	///
	/// \brief	A cube centered in (0,0,0), 1.0 height/width....
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Cube : public Geometry
	{
	protected:
		Square * m_square[6] ;

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Cube::Cube(Material * material)
		///
		/// \brief	Constructeur
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	material	The material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Cube(Material * material)
			: Geometry()
		{
			Square sq0(material) ;
			sq0.translate(Math::makeVector(0.0f,0.0f,0.5f)) ;
			merge(sq0) ;

			Square sq1(material) ;
			sq1.translate(Math::makeVector(0.0f,0.0f,0.5f)) ;
			sq1.rotate(Math::Quaternion<double>(Math::makeVector(1.0f,0.0f,0.0f), (double)M_PI/2.0f)) ;
			merge(sq1) ;

			Square sq2(material) ;
			sq2.translate(Math::makeVector(0.0f,0.0f,0.5f)) ;
			sq2.rotate(Math::Quaternion<double>(Math::makeVector(1.0f,0.0f,0.0f), (double)M_PI)) ;
			merge(sq2) ;
			
			Square sq3(material) ;
			sq3.translate(Math::makeVector(0.0f,0.0f,0.5f)) ;
			sq3.rotate(Math::Quaternion<double>(Math::makeVector(1.0f,0.0f,0.0f), (double)-M_PI/2.0f)) ;
			merge(sq3) ;

			Square sq4(material) ;
			sq4.translate(Math::makeVector(0.0,0.0,0.5)) ;
			sq4.rotate(Math::Quaternion<double>(Math::makeVector(0.0f,1.0f,0.0f), (double)M_PI/2.0f)) ;
			merge(sq4) ;

			Square sq5(material) ;
			sq5.translate(Math::makeVector(0.0f,0.0f,0.5f)) ;
			sq5.rotate(Math::Quaternion<double>(Math::makeVector(0.0f,1.0f,0.0f), (double)-M_PI/2.0f)) ;
			merge(sq5) ;

			// Déclaration des carrés constituant les faces
			//m_square[0] = new Square(material) ;
			//m_square[0]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[1] = new Square(material) ;
			//m_square[1]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[1]->rotate(Math::Quaternion(Math::Vector3(1.0f,0.0f,0.0f), (double)M_PI/2.0f)) ;
			//m_square[2] = new Square(material) ;
			//m_square[2]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[2]->rotate(Math::Quaternion(Math::Vector3(1.0f,0.0f,0.0f), (double)M_PI)) ;
			//m_square[3] = new Square(material) ;
			//m_square[3]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[3]->rotate(Math::Quaternion(Math::Vector3(1.0f,0.0f,0.0f), (double)-M_PI/2.0f)) ;
			//m_square[4] = new Square(material) ;
			//m_square[4]->translate(Math::Vector3(0.0,0.0,0.5)) ;
			//m_square[4]->rotate(Math::Quaternion(Math::Vector3(0.0f,1.0f,0.0f), (double)M_PI/2.0f)) ;
			//m_square[5] = new Square(material) ;
			//m_square[5]->translate(Math::Vector3(0.0f,0.0f,0.5f)) ;
			//m_square[5]->rotate(Math::Quaternion(Math::Vector3(0.0f,1.0f,0.0f), (double)-M_PI/2.0f)) ;
			// Ajout des géométries dans cette géométrie
			//merge(m_square[0]) ;
			//merge(m_square[1]) ;
			//merge(m_square[2]) ;
			//merge(m_square[3]) ;
			//merge(m_square[4]) ;
			//merge(m_square[5]) ;
		}
	} ;
}

#endif
