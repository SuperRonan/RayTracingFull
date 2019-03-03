#ifndef _Geometry_Square_H
#define _Geometry_Square_H

#include <Geometry/Geometry.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Square
	///
	/// \brief	A square on plane (X,Y), centered in (0,0,0).
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Square : public Geometry
	{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Square::Square(Material * material)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	material	If non-null, the material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Square(Material * material)
			: Geometry()
		{
			int p0 = addVertex(Math::makeVector(0.5, 0.5, 0.0)) ;
			int p1 = addVertex(Math::makeVector(0.5, -0.5, 0.0)) ;
			int p2 = addVertex(Math::makeVector(-0.5, 0.5, 0.0)) ;
			int p3 = addVertex(Math::makeVector(-0.5, -0.5, 0.0)) ;
			addTriangle(p0,p1,p2,material) ;
			addTriangle(p3,p1,p2, material) ;
		}
	} ;





	class OneTriangle : public Geometry
	{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Square::Square(Material * material)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	material	If non-null, the material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		OneTriangle(Material * material)
			: Geometry()
		{
			int p0 = addVertex(Math::makeVector(0.5, 0.5, 0.0));
			int p1 = addVertex(Math::makeVector(0.5, -0.5, 0.0));
			int p2 = addVertex(Math::makeVector(-0.5, 0.5, 0.0));
			addTriangle(p0, p1, p2, material);
		}
	};
}

#endif
