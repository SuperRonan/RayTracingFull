#ifndef _Geometry_camera_H
#define _Geometry_camera_H

#include <Math/Vectorf.h>
#include <Math/Quaternion.h>
#include <Geometry/Ray.h>
#include <math.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Camera
	///
	/// \brief	A camera.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Camera
	{
	public:
		/// \brief	The camera position.
		Math::Vector3f m_position ;
		/// \brief	The aim of the camera.
		Math::Vector3f m_target ;
		/// \brief	Distance of the focal plane.
		double		  m_planeDistance ;
		/// \brief	Width of the projection rectangle.
		double	      m_planeWidth ;
		/// \brief	Height of the projection rectangle.
		double		  m_planeHeight ;
		/// \brief	The front vector of the camera.
		Math::Vector3f m_front ;
		/// \brief	The right vector.
		Math::Vector3f m_right ;
		/// \brief	The down vector.
		Math::Vector3f m_down ;
		/// \brief	The width vector of the projection rectangle.
		Math::Vector3f m_widthVector ;
		/// \brief	The height vector of the projection rectangle.
		Math::Vector3f m_heightVector ;
		/// \brief	The upper left point oft he projection rectangle.
		Math::Vector3f m_upLeftPoint ;
	
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Camera::computeParameters()
		///
		/// \brief	Calculates the camera parameters.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void computeParameters()
		{
			m_front = m_target-m_position ;
			m_front = m_front/m_front.norm() ;
			m_right = Math::Quaternion<double>(Math::makeVector(0.0f, 0.0f, 1.0f), -3.14159265f/2.0f).rotate(m_front) ;
			m_right[2] = 0.0f;
			m_right = m_right/m_right.norm() ;
			m_down  = m_front ^ m_right ;
			m_down  = m_down/m_down.norm() ;
			m_widthVector  = m_right*m_planeWidth ;
			m_heightVector = m_down*m_planeHeight ;
			m_upLeftPoint  = m_position+m_front*m_planeDistance-m_widthVector*0.5-m_heightVector*0.5 ;
		}

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Camera::Camera(Math::Vector3 const & position = Math::Vector3(0.0, 0.0, 0.0),
		/// 	Math::Vector3 const & target = Math::Vector3(0.0, 1.0, 0.0), double planeDistance=1.0,
		/// 	double planeWidth=1.0, double planeHeight=1.0)
		///
		/// \brief	Constructeur de la caméra.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	position	 	The camera position
		/// \param	target		 	The target of the camera
		/// \param	planeDistance	La distance of the focal plane.
		/// \param	planeWidth   	Width of the projection rectangle.
		/// \param	planeHeight  	Height of the projection rectangle.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Camera(Math::Vector3f const & position = Math::makeVector(0.0f, 0.0f, 0.0f), 
			   Math::Vector3f const & target = Math::makeVector(0.0f, 1.0f, 0.0f), 
			   double planeDistance=1.0f, double planeWidth=1.0f, double planeHeight=1.0f)
		   : m_position(position), m_target(target), m_planeDistance(planeDistance), 
		     m_planeWidth(planeWidth), m_planeHeight(planeHeight)
		{
			computeParameters() ;
		}

		/// <summary>
		/// Translates the camera in local coordinates (X = right, Y = front, Z=up).
		/// </summary>
		/// <param name="translation">The translation vector.</param>
		void translateLocal(Math::Vector3f const & translation)
		{
			Math::Vector3f trans =m_right*translation[0] + m_front*translation[1] - m_down*translation[2];
			m_position = m_position + trans;
			m_target = m_target + trans;
			computeParameters();
		}


		void update_both(Math::Vector3f const & position, Math::Vector3f const & target)
		{
			m_position = position;
			m_target = target;
			computeParameters();
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Camera::setPosition(Math::Vector3 const & position)
		///
		/// \brief	Sets the camera position position.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	position	The new camera position.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setPosition(Math::Vector3f const & position)
		{
			m_position = position ;
			computeParameters() ;
		}

		Math::Vector3f const& getPosition()
		{
			return m_position;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Camera::setTarget(Math::Vector3 const & target)
		///
		/// \brief	Sets the target.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	target	The new target.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setTarget(Math::Vector3f const & target)
		{
			m_target = target ;
			computeParameters() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Ray Camera::getRay(double coordX, double coordY) const
		///
		/// \brief	Get a primary ray from screen coordinates (coordX, coordY)
		/// 		
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	coordX	X coordinate in the projection rectangle.
		/// \param	coordY	Y coordinate in the projection rectangle.
		///
		/// \return	The ray.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Ray getRay(double coordX, double coordY) const
		{
			return Ray(m_position, m_upLeftPoint+m_widthVector*coordX+m_heightVector*coordY-m_position) ;
		}
	} ;
}

#endif
