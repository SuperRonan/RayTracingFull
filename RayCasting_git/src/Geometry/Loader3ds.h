#ifndef _HelperGl_Loader3ds_H
#define _HelperGl_Loader3ds_H

#include <lib3ds/file.h>
#include <lib3ds/mesh.h>
#include <lib3ds/material.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <Geometry/Material.h>
#include <Geometry/Geometry.h>

namespace Geometry
{
	class Loader3ds
	{
	private:
		/// \brief The 3DS file name
		Lib3dsFile * m_file;
		/// \brief The loaded materials
		::std::map<::std::string, Material *> m_materials;
		/// \brief The loaded meshes
		::std::vector<Geometry*> m_meshes ;


		/// <summary>
		/// Loads a material.
		/// </summary>
		/// <param name="material">The material.</param>
		/// <param name="texturePath">The texture path.</param>
		void loadMaterial(Lib3dsMaterial * material, const ::std::string & texturePath);

		/// <summary>
		/// Loads the materials.
		/// </summary>
		/// <param name="texturePath">The texture path.</param>
		void loadMaterials(const ::std::string & texturePath);

		/// <summary>
		/// Loads a mesh.
		/// </summary>
		/// <param name="mesh">The mesh.</param>
		void loadMesh(Lib3dsMesh * mesh);

		/// <summary>
		/// Loads the meshes.
		/// </summary>
		void loadMeshes()
		{
			for(Lib3dsMesh * mesh = m_file->meshes ; mesh!=NULL ; mesh = mesh->next)
			{
				loadMesh(mesh) ;
			}
		}

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Loader3ds::Loader3ds(const ::std::string & filename, const ::std::string & texturePath);
		///
		/// \brief	Constructor which loads the provided 3ds file.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	08/02/2016
		///
		/// \param	filename   	Filename of the file.
		/// \param	texturePath	Full pathname of the texture file.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Loader3ds(const ::std::string & filename, const ::std::string & texturePath);

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const ::std::vector<Geometry*> & Loader3ds::getMeshes()
		///
		/// \brief	Gets the loaded meshes.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \return	The meshes.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const ::std::vector<Geometry*> & getMeshes() const;

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	::std::vector<Material*> Loader3ds::getMaterials()
		///
		/// \brief	Gets the loaded materials.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		////////////////////////////////////////////////////////////////////////////////////////////////////
		::std::vector<Material*> getMaterials() const;
	};
}

#endif 