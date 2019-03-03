#ifndef _Geometry_RayTriangleIntersection_H
#define _Geometry_RayTriangleIntersection_H

#include <Geometry/Ray.h>
#include <Geometry/Triangle.h>
#include <Spy/Spy.h>
#include <assert.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	RayTriangleIntersection
	///
	/// \brief	An intersection between a ray and a triangle.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class RayTriangleIntersection
	{
	protected:
		/// \brief	The distance between ray source and the intersection.
		double m_t ;
		/// \brief	The u coordinate of the intersection.
		double m_u ;
		/// \brief	The v coordinate of the intersection.
		double m_v ;
		/// \brief	Is the intersection valid?
		bool m_valid ;
		/// \brief	The triangle associated to the intersection.
		const Triangle * m_triangle ;

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RayTriangleIntersection::RayTriangleIntersection(const Triangle * triangle,
		/// 	const Ray & ray)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	triangle	The triangle.
		/// \param	ray			The ray.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RayTriangleIntersection(const Triangle * triangle, const Ray & ray)
			: m_triangle(triangle)
		{
			m_valid=triangle->intersection(ray, m_t, m_u, m_v) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RayTriangleIntersection::RayTriangleIntersection()
		///
		/// \brief	Constructor that constructs an invalid intersection.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RayTriangleIntersection()
			: m_valid(false), m_triangle(NULL)
		{}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool RayTriangleIntersection::valid() const
		///
		/// \brief	Is the intersection valid?
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	true if valid, false otherwise.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool valid() const
		{ return m_valid ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	double RayTriangleIntersection::tRayValue() const
		///
		/// \brief	Distance between the ray source and the intersection point
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		double tRayValue() const
		{ return m_t ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	double RayTriangleIntersection::uTriangleValue() const
		///
		/// \brief	The u value of the intersection.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		double uTriangleValue() const
		{ return m_u ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	double RayTriangleIntersection::vTriangleValue() const
		///
		/// \brief	The v value of the intersection
		/// 		
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		double vTriangleValue() const
		{ return m_v ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Triangle * RayTriangleIntersection::triangle() const
		///
		/// \brief	Returns the triangle associated with the intersection.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	The triangle.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Triangle * triangle() const
		{ return m_triangle ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 RayTriangleIntersection::intersection() const
		///
		/// \brief	Returns the intersection point.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector3f intersection() const
		{
			return m_triangle->samplePoint(m_u, m_v);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool RayTriangleIntersection::operator< (RayTriangleIntersection const & i) const
		///
		/// \brief	Compares two intersections associated with the same ray (otherwise, it is a non-sense)
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	i	Zero-based index of the.
		///
		/// \return	true if the first parameter is less than the second.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool operator < (RayTriangleIntersection const & i) const
		{ 
			return (m_valid & i.m_valid & (m_t<i.m_t)) | (!i.m_valid) ; 
		}
	} ;
}

#endif
