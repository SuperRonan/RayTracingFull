#include <Geometry/Loader3ds.h>
#include <lib3ds/mesh.h>

namespace Geometry
{
	void Loader3ds::loadMaterial( Lib3dsMaterial * material, const ::std::string & texturePath )
	{
		Material currentMaterial ;
		RGBColor color ;
		color.set(material->diffuse);
		currentMaterial.setDiffuse(color);
		color.set(material->specular);
		currentMaterial.setSpecular(color);
		color.set(material->ambient);
		currentMaterial.setAmbient(color);
		currentMaterial.setShininess(material->shininess);
		::std::string textureName = material->texture1_map.name;
		if(textureName!="")
		{
			currentMaterial.setTextureFile(texturePath+"\\"+textureName);
		}
		::std::string materialName = material->name;
		m_materials[materialName] = new Material(currentMaterial);
		::std::cout<<"Loader3ds: created material "<<materialName<<" "<<currentMaterial.getDiffuse()<<currentMaterial.getSpecular()<<", texture: "<< textureName <<::std::endl ;
	}

	void Loader3ds::loadMaterials( const ::std::string & texturePath )
	{
		for(Lib3dsMaterial * material = m_file->materials ; material != NULL ; material = material->next)
		{
			loadMaterial(material, texturePath);
		}
	}

	void Loader3ds::loadMesh( Lib3dsMesh * mesh )
	{
		//::std::string meshName = mesh->name;
		::std::cout << "Loading mesh " << mesh->name << ::std::flush;
		Geometry * geometry = new Geometry;
		// 1.1 - We register the vertices
		for (unsigned int cpt = 0; cpt < mesh->points; ++cpt)
		{
			Lib3dsPoint pt = mesh->pointL[cpt];
			Math::Vector3f v = Math::makeVector(pt.pos[0], pt.pos[1], pt.pos[2]);
			geometry->addVertex(v);
			//vertices.push_back(v);
		}
		// 1.2 - We register the texture coordinates (if any)
		if (mesh->texels > 0)
		{
			::std::cout << ", has texture coordinates" << ::std::flush;
		}
		for(unsigned int cpt=0 ; cpt<mesh->texels ; ++cpt)
		{
			Math::Vector2f t = Math::makeVector(mesh->texelL[cpt][0], mesh->texelL[cpt][1]);
			geometry->addTextureCoordinates(t);
		}
		// 1.3 we compute the normals
		Lib3dsVector * normals = new Lib3dsVector[3*mesh->faces];
		Math::Vector3f * vertexNormals = new Math::Vector3f[3 * mesh->faces];
		lib3ds_mesh_calculate_normals(mesh, normals);
		for (unsigned int cpt = 0; cpt < 3 * mesh->faces; ++cpt)
		{
			vertexNormals[cpt] = Math::makeVector(normals[cpt][0], normals[cpt][1], normals[cpt][2]).normalized();
		}
		// 2 - We register the triangles
		for (unsigned int cpt = 0; cpt < mesh->faces; ++cpt)
		{
			Lib3dsFace tmp = mesh->faceL[cpt];
			if (m_materials.find(tmp.material) == m_materials.end())
			{
				::std::cout << "Problem, material " << tmp.material << " not found..." << ::std::endl;
			}
			else
			{
				geometry->addTriangle(tmp.points[0], tmp.points[1], tmp.points[2], m_materials[tmp.material], vertexNormals+(cpt*3));
			}
		}
		// 3 - We register the mesh
		m_meshes.push_back(geometry);
		// 4 - Cleanup
		delete[] normals;
		::std::cout << std::endl;
	}

	Loader3ds::Loader3ds( const ::std::string & filename, const ::std::string & texturePath )
	{
		::std::cout<<"Loader3ds: loading file "<<filename<<::std::endl; 
		m_file = lib3ds_file_load(filename.c_str()) ;
		if(m_file==NULL)
		{
			::std::cerr<<"Loader3ds: unable to load file "+filename<<::std::endl ;
			return ;
		}
		loadMaterials(texturePath) ;
		loadMeshes();
		::std::cout<<"Loader3ds: OK"<<::std::endl ;
	}

	const ::std::vector<Geometry*> & Loader3ds::getMeshes() const
	{
		return m_meshes ;
	}

	::std::vector<Material*> Loader3ds::getMaterials() const
	{
		::std::vector<Material*> result ;
		::std::transform(m_materials.begin(), m_materials.end(), ::std::back_inserter(result), [](::std::pair<::std::string, Material*> const & m) -> Material* { return m.second ; } );
		return result;
	}

}