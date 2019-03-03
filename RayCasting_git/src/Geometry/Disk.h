#ifndef _Geometry_Disk_H
#define _Geometry_Disk_H

#include <Geometry/Geometry.h>
#include <Geometry/Material.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Disk
	///
	/// \brief	A disk on (X,Y) plane, centered in (0,0,0) with a 0.5 radius.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	///////////////////////////////////////////////////////////////////////////////////////////////////
	class Disk : public Geometry
	{
	protected:
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Disk::Disk(int nbDiv, Material * material)
		///
		/// \brief	Construtor
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	nbDiv   	Number of subdivisions of the circle delimiting the disk.
		/// \param	material	The material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Disk(int nbDiv, Material * material)
		{
			unsigned int center = addVertex(Math::Vector3f()) ;
			::std::vector<unsigned int> vertices ;
			for(int cpt=0 ; cpt<nbDiv ; cpt++)
			{
				double angle = double((2.0f*M_PI/nbDiv)*cpt) ;
				int i = addVertex(Math::makeVector(cos(angle), sin(angle), 0.0)) ;
				vertices.push_back(i) ;
			}
			for(int cpt=0 ; cpt<nbDiv ; cpt++)
			{
				addTriangle(center, vertices[cpt], vertices[(cpt+1)%nbDiv], material) ;
			}

			//Math::Vector3 * center = new Math::Vector3();
			//// Generation des points
			//for(int cpt=0 ; cpt<nbDiv ; cpt++)
			//{
			//	double angle = double((2.0f*M_PI/nbDiv)*cpt) ;
			//	addVertex(new Math::Vector3(cos(angle), sin(angle), 0.0)) ;
			//}
			//addVertex(center) ;
			// Generation des triangles
			//for(int cpt=0 ; cpt<nbDiv ; cpt++)
			//{
			//	Math::Vector3 * pt1 = m_vertices[(cpt)%(nbDiv)] ;
			//	Math::Vector3 * pt2 = m_vertices[(cpt+1)%(nbDiv)] ;
			//	addTriangle(new Triangle(center, pt1, pt2, material)) ;
			//}
		}
	} ;
} ;

#endif
