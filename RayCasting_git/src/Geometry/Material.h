#ifndef _Geometry_Material_H
#define _Geometry_Material_H

#include <Geometry/RGBColor.h>
#include <Geometry/Texture.h>
#include "medium.h"

namespace Geometry
{
	/** \brief A material definition. */
	class Material
	{
	protected:
		/// <summary> The ambient color </summary>
		RGBColor m_ambientColor ;
		/// <summary> the diffuse color</summary>
		RGBColor m_diffuseColor ;
		/// <summary> The specular color</summary>
		RGBColor m_specularColor ;
		/// <summary> The shininess</summary>
		double    m_shininess ;
		/// <summary> The emissive color</summary>
		RGBColor m_emissiveColor ;
		/// <summary> The filename of the itexture image</summary>
		std::string m_textureFile;
		/// <summary> The texture or nullptr if no texture is associated with this material</summary>
		Texture * m_texture; 

		

	public:
		const medium * m_medium = nullptr;

		bool m_transparent = false;

		/// <summary>
		/// Initializes a new instance of the <see cref="Material"/> class.
		/// </summary>
		/// <param name="ambientColor">The ambient color.</param>
		/// <param name="diffuseColor">The diffuse color.</param>
		/// <param name="specularColor">The specular color.</param>
		/// <param name="shininess">The shininess.</param>
		/// <param name="emissiveColor">The emissive color.</param>
		Material(RGBColor const & ambientColor = RGBColor(), RGBColor const & diffuseColor = RGBColor(),
				 RGBColor specularColor = RGBColor(), double shininess = 1.0, RGBColor const & emissiveColor = RGBColor(), 
				 std::string const & textureFile = "")
				 : m_ambientColor(ambientColor), m_diffuseColor(diffuseColor), m_specularColor(specularColor),
				   m_shininess(shininess), m_emissiveColor(emissiveColor), m_textureFile(textureFile), m_texture(NULL)
		{}

		/// <summary>
		/// Sets the ambient color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setAmbient(RGBColor const & color)
		{
			m_ambientColor = color;
		}

		/// <summary>
		/// Gets the ambient color.
		/// </summary>
		/// <returns>The ambiant color</returns>
		const RGBColor & getAmbient() const
		{ return m_ambientColor ; }

		/// <summary>
		/// Sets the diffuse color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setDiffuse(RGBColor const & color)
		{
			m_diffuseColor = color;
		}
		
		/// <summary>
		/// Gets the diffuse color.
		/// </summary>
		/// <returns></returns>
		const RGBColor & getDiffuse() const
		{ return m_diffuseColor ; }

		/// <summary>
		/// Sets the specular color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setSpecular(RGBColor const & color)
		{
			m_specularColor = color;
		}

		/// <summary>
		/// Gets the specular color.
		/// </summary>
		/// <returns></returns>
		const RGBColor & getSpecular() const
		{ return m_specularColor ; }

		/// <summary>
		/// Sets the shininess.
		/// </summary>
		/// <param name="s">The shininess.</param>
		void setShininess(double s)
		{
			m_shininess = s;
		}

		/// <summary>
		/// Gets the shininess.
		/// </summary>
		/// <returns></returns>
		const double & getShininess() const
		{ return m_shininess ; }

		/// <summary>
		/// Sets the emissive color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setEmissive(RGBColor const & color)
		{
			m_emissiveColor = color;
		}

		/// <summary>
		/// Gets the emissive color.
		/// </summary>
		/// <returns></returns>
		const RGBColor & getEmissive() const
		{ return m_emissiveColor ; }

		/// <summary>
		/// Sets the texture file.
		/// </summary>
		/// <param name="textureFile">The texture file.</param>
		void setTextureFile(const ::std::string & textureFile)
		{
			m_textureFile = textureFile;
			::std::cout << "Loading texture: "<< m_textureFile << "..." << ::std::flush;
			m_texture = new Texture(m_textureFile);
			if (!m_texture->isValid())
			{
				delete m_texture;
				m_texture = NULL;
				::std::cout << "discarded" << ::std::endl;
			}
			else
			{
				::std::cout << "OK" << ::std::endl;
			}
		}

		/// <summary>
		/// Gets the texture file.
		/// </summary>
		/// <returns></returns>
		const ::std::string & getTextureFile() const
		{
			return m_textureFile;
		}

		/// <summary>
		/// Returns the texture
		/// </summary>
		/// <returns></returns>
		const Texture & getTexture() const
		{
			return *m_texture;
		}

		/// <summary>
		/// Tests if a texture is associated with this material.
		/// </summary>
		/// <returns> true if this materail has a texture, false otherwise</returns>
		bool hasTexture() const
		{
			return m_texture != NULL;
		}




		bool is_transparent()const
		{
			return m_transparent;
		}
	};
}

#endif
