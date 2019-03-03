#ifndef _Geometry_BoundingBox_H
#define _Geometry_BoundingBox_H

#include <limits>
#include <Geometry/Geometry.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	BoundingBox
	///
	/// \brief	A bounding box.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class BoundingBox
	{
	protected:
		/// \brief	The bounds (min / max vectors).
		Math::Vector3f m_bounds[2] ;
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	BoundingBox::BoundingBox(Geometry const & geometry)
		///
		/// \brief	Initializes the bounding box from the provided geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	geometry	The geometry.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		BoundingBox(Geometry const & geometry)
		{
			set(geometry);
		}




#define min_val std::numeric_limits<double>::lowest()
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	BoundingBox::BoundingBox(Math::Vector3 const & minVertex, Math::Vector3 const & maxVertex)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	minVertex	The smallest coordinates on X, Y, Z axes.
		/// \param	maxVertex	The highest coordinates on X, Y, Z axes.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		BoundingBox(Math::Vector3f const & minVertex = Math::makeVector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()), 
					Math::Vector3f const & maxVertex = Math::makeVector(min_val, min_val, min_val))
		{
			m_bounds[0] = minVertex ;
			m_bounds[1] = maxVertex ;
		}

		/// <summary>
		/// Tests wether the box is empty of not.
		/// </summary>
		/// <returns></returns>
		bool isEmpty() const
		{
			bool result = false;
			for (int cpt = 0; cpt < 3; ++cpt)
			{
				result |= m_bounds[0][cpt] > m_bounds[1][cpt];
			}
			return result;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void BoundingBox::set( Geometry const &geometry )
		///
		/// \brief	Sets the bouding box with the given geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	10/12/2013
		///
		/// \param	geometry	The geometry.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void set( Geometry const &geometry ) 
		{
			const std::deque<Math::Vector3f> & vertices(geometry.getVertices()) ;
			m_bounds[0] = vertices[0] ;
			m_bounds[1] = vertices[0] ;
			update(geometry) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void BoundingBox::update(Geometry const & geometry)
		///
		/// \brief	Updates the bounding box with the given geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	10/12/2013
		///
		/// \param	geometry	The geometry.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void update(Geometry const & geometry)
		{
			const std::deque<Math::Vector3f> & vertices(geometry.getVertices()) ;
			for(auto it=vertices.begin(), end=vertices.end() ; it!=end ; ++it)
			{
				m_bounds[0] = m_bounds[0].simdMin(*it) ;
				m_bounds[1] = m_bounds[1].simdMax(*it) ;
			}
		}

		/// <summary>
		/// Updates the bounding box with the provided triangle.
		/// </summary>
		/// <param name="triangle"></param>
		void update(const Triangle & triangle)
		{
			update(triangle.vertex(0));
			update(triangle.vertex(1));
			update(triangle.vertex(2));
		}

		/// <summary>
		/// Updates the bounding box with the provided point.
		/// </summary>
		/// <param name="v"></param>
		void update(const Math::Vector3f & v)
		{
			m_bounds[0] = m_bounds[0].simdMin(v);
			m_bounds[1] = m_bounds[1].simdMax(v);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void BoundingBox::update(BoundingBox const & boundingBox)
		///
		/// \brief	Updates this bounding box to bound the given boundingBox.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	10/12/2013
		///
		/// \param	boundingBox	The bounding box.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void update(BoundingBox const & boundingBox)
		{
			m_bounds[0] = m_bounds[0].simdMin(boundingBox.m_bounds[0]) ;
			m_bounds[1] = m_bounds[1].simdMax(boundingBox.m_bounds[1]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool BoundingBox::intersect(const Ray & ray) const
		///
		/// \brief	Tests if the provided ray intersects this box.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	09/12/2013
		///
		/// \param	ray	The ray.
		///
		/// \return	true if an intersection is found, false otherwise.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool intersect(const Ray & ray, double t0, double t1, double & entryT, double & exitT) const
		{
			//int sign[3] = { ray.direction()[0]<0.0, ray.direction()[1]<0.0, ray.direction()[2]<0.0 } ;
			const int * sign = ray.getSign() ;
			Math::Vector3f tmin = Math::makeVector(m_bounds[sign[0]][0], m_bounds[sign[1]][1], m_bounds[sign[2]][2]) ;
			tmin = (tmin - ray.source()).simdMul(ray.invDirection()) ;
			Math::Vector3f tmax = Math::makeVector(m_bounds[1-sign[0]][0], m_bounds[1-sign[1]][1], m_bounds[1-sign[2]][2]) ;
			tmax = (tmax - ray.source()).simdMul(ray.invDirection()) ;
			
			if((tmin[0]>tmax[1]) || (tmin[1]>tmax[0]))
			{
				return false ;
			}
			if(tmin[1] > tmin[0])
			{
				tmin[0] = tmin[1] ;
			}
			if(tmax[1] < tmax[0])
			{
				tmax[0] = tmax[1] ;
			}
			
			if((tmin[0]>tmax[2]) || (tmin[2] > tmax[0]))
			{
				return false ;
			}
			if(tmin[2] > tmin[0])
			{
				tmin[0] = tmin[2] ;
			}
			if(tmax[2] < tmax[0])
			{
				tmax[0] = tmax[2] ; 
			}
			bool intersectionFound = (tmin[0]<t1) && (tmax[0]>t0) ;
			if (intersectionFound)
			{
				entryT = ::std::max(tmin[0],t0);
				exitT = ::std::min(tmax[0],t1);
			}
			return intersectionFound;
		}

		/// <summary>
		/// The min vector.
		/// </summary>
		/// <returns></returns>
		Math::Vector3f min() const
		{
			return m_bounds[0];
		}

		Math::Vector3f& min()
		{
			return m_bounds[0];
		}

		/// <summary>
		/// The max vector.
		/// </summary>
		/// <returns></returns>
		Math::Vector3f max() const
		{
			return m_bounds[1];
		}

		Math::Vector3f & max() 
		{
			return m_bounds[1];
		}

		double surface()const
		{
			Math::Vector3f diag = max() - min();
			return 2.0 * (diag * diag);
		}

		bool operator==(BoundingBox const & other)
		{
			return min() == other.min() && max() == other.max();
		}

		bool operator!=(BoundingBox const & other)
		{
			return min() != other.min() || max() != other.max();
		}
	} ;
}

#endif