#ifndef _Geometry_Scene_H
#define _Geometry_Scene_H

#include <windows.h>
#include <Geometry/Geometry.h>
#include <Geometry/PointLight.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Camera.h>
#include <Geometry/BoundingBox.h>
#include <Math/RandomDirection.h>
#include <windows.h>
#include <System/aligned_allocator.h>
#include <Math/Constant.h>
#include <queue>
#include <functional>
#include <random>
#include <Geometry/LightSampler.h>
#include <cmath>
#include <limits>
#include "medium.h"
#include <stack>
#include <utils.h>
#include <valarray>
#include <tbb/parallel_for_each.h>
#include <functional>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Scene
	///
	/// \brief	An instance of a geometric scene that can be rendered using ray casting. A set of methods
	/// 		allowing to add geometry, lights and a camera are provided. Scene rendering is achieved by
	/// 		calling the Scene::compute method.
	///			A scene doesnt implement any acceleration stuctures, exept for testing the intersection between
	///			the ray and a geomtry's bounding box
	///
	/// \author	F. Lamarche, Université de Rennes 1 and Ronan Cailleau
	/// \date	03/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Scene
	{
	public:
		/// \brief	The visualizer (rendering target).
		Visualizer::Visualizer * m_visu ;
		/// \brief	The scene geometry (basic representation without any optimization).
		::std::deque<::std::pair<BoundingBox, Geometry> > m_geometries;
		//Geometry m_geometry ;
		/// \brief	The lights.
		std::vector<PointLight> m_lights;
		/// \brief	The camera.
		Camera m_camera;
		/// \brief The scene bounding box
		BoundingBox m_sceneBoundingBox;
		/// \brief number of diffuse samples for global illumination
		size_t m_diffuseSamples;
		/// \brief Number of specular samples 
		size_t m_specularSamples;
		/// \brief Number of light samples 
		size_t m_lightSamples;
		/// brief Rendering pass number
		int m_pass;
		/// <summary>
		/// The light sampler associated with the scene
		/// </summary>
		std::vector<LightSampler> m_lightSamplers;

		RGBColor m_background_color;

		Texture m_skybox;

		medium m_medium;

		std::function<RGBColor (Scene const*, RayTriangleIntersection const&, Ray const&, Math::Vector3f &, Math::Vector3f &, Math::Vector3f &)> m_phong;


		size_t m_number_triangle = 0;

		

	public:

		size_t max_depth = 5;


		


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Scene::Scene(Visualizer::Visualizer * visu)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \param [in,out]	visu	If non-null, the visu.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Scene(Visualizer::Visualizer * visu, std::string const& skybox_file="..\\..\\skybox.png", medium med=medium())
			: m_visu(visu), 
			m_diffuseSamples(30), 
			m_specularSamples(30),
			m_lightSamples(0), 
			m_background_color(0.1, 0.1, 0.1),
			m_skybox(skybox_file), 
			m_medium(med),
			m_phong(&Scene::phong_ray)
		{}

		/// <summary>
		/// Prints stats about the geometry associated with the scene
		/// </summary>
		void printStats()
		{
			::std::cout << "Scene: " << m_number_triangle << " triangles" << ::std::endl;
		}

		/// <summary>
		/// Computes the scene bounding box.
		/// </summary>
		/// <returns></returns>
		const BoundingBox & getBoundingBox()
		{
			return m_sceneBoundingBox;
		}

		/// <summary>
		/// Sets the number of diffuse samples
		/// </summary>
		/// <param name="number"> The number of diffuse samples</param>
		void setDiffuseSamples(size_t number)
		{
			m_diffuseSamples = number;
		}

		/// <summary>
		/// Sets the number of specular samples
		/// </summary>
		/// <param name="number"></param>
		void setSpecularSamples(size_t number)
		{
			m_specularSamples = number;
		}

		/// <summary>
		/// Sets the number of light samples if the scene contains surface area lights
		/// </summary>
		/// <param name="number">The number of samples</param>
		void setLightSamples(size_t number)
		{
			m_lightSamples = number;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::add(const Geometry & geometry)
		///
		/// \brief	Adds a geometry to the scene.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \param	geometry The geometry to add.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void add(const Geometry & geometry)
		{
			if (geometry.getVertices().size() == 0) { return; }
			BoundingBox box(geometry) ;
			m_geometries.push_back(::std::make_pair(box, geometry)) ;
			m_geometries.back().second.computeVertexNormals(Math::piDiv4/2);
			if (m_geometries.size() == 1)
			{
				m_sceneBoundingBox = box;
			}
			else
			{
				m_sceneBoundingBox.update(box);
			}
			m_number_triangle += geometry.getTriangles().size();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::add(PointLight * light)
		///
		/// \brief	Adds a poitn light in the scene.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	light	If non-null, the light to add.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void add(const PointLight & light)
		{
			m_lights.push_back(light) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::setCamera(Camera const & cam)
		///
		/// \brief	Sets the camera.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	cam	The camera.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setCamera(Camera const & cam)
		{
			m_camera = cam ;
		}

		/////////////////////////////////////////////////////////////
		// Computes the intersecion between the ray and the scene
		// intersection with the scene, exept the current triangle
		/////////////////////////////////////////////////////////////
		virtual bool intersection(Ray const& ray, RayTriangleIntersection & intersect_info, const Triangle * current=nullptr)const
		{
			CastedRay cray(ray);
			double current_min_t = 0, current_max_t = std::numeric_limits<double>::max();
			for (const std::pair< BoundingBox, Geometry > & geom : m_geometries)
			{
				double t1, t2;
				if (geom.first.intersect(cray, current_min_t, current_max_t, t1, t2))
				{
					//current_min_t = ::std::min(t1, std::min(t2, current_min_t));
					
					//geom.second.intersection(cray);
					for (Triangle const& t : geom.second.getTriangles())
					{
						if (&t != current)
						{
							cray.intersect(&t);
						}
					}
					if (cray.validIntersectionFound())
					{
						
						current_max_t = std::min(current_max_t, cray.intersectionFound().tRayValue());
					}
					
				}


				
			}
			if (cray.validIntersectionFound())
			{
				intersect_info = cray.intersectionFound();
				return true;
			}
			else
			{
				return false;
			}
		}

	protected:


		



	public:


		/////////////////////////////////////////////////////
		// Computes the phong: with sending a ray to each light
		/////////////////////////////////////////////////////
		RGBColor phong_ray(RayTriangleIntersection const& rti, Ray const& ray, Math::Vector3f & normal, Math::Vector3f & to_view, Math::Vector3f & reflected)const
		{
			RGBColor res(0.0, 0.0, 0.0);
			Math::Vector3f const& intersection = rti.intersection();
			Triangle const& tri = *(rti.triangle());
			Material & mat = *(tri.material());

			double u = rti.uTriangleValue(), v = rti.vTriangleValue();

			RGBColor const tex_color = tri.sampleTexture(u, v);

			double dist_inter_source = rti.tRayValue();;

			res = res + mat.getAmbient();
			res = res + mat.getEmissive();

			//res = res / (dist_inter_source + 1);

			to_view = (ray.source() - intersection).normalized();

			normal = tri.sampleNormal(u, v, ray.source());

			

			RGBColor mat_diffuse = mat.getDiffuse();
			RGBColor mat_specular = mat.getSpecular();

			reflected = normal * (2 * (normal * to_view)) - to_view;

			//TriangleMultipleTexture<bool> const& light_texture = tri.m_light_tex;

			if (!(mat_diffuse.isBlack() && mat_specular.isBlack()))
			{

				for (unsigned int i=0; i<m_lights.size(); ++i)
				{
					PointLight const& l = m_lights[i];
					Math::Vector3f to_light = (l.position() - intersection).normalized();

					if (to_light * normal < 0.f)
					{
						continue;
					}

					//RGBColor incoming_light = send_ray_to_light(rti, l); // light to tri
					
					
					RGBColor incoming_light = send_ray_to_light(Ray(intersection, to_light), l);//tri to light
				
					

					if (incoming_light.isBlack())
					{
						continue;
					}
					
					//Diffuse
					if (!mat_diffuse.isBlack())
					{
						RGBColor diffuse = mat_diffuse * (to_light * normal);
						res = res + diffuse * incoming_light;
					}
					
					//Specular
					if (!mat_specular.isBlack())
					{
						
						if (reflected * to_light >= 0)
						{


							RGBColor specular = mat_specular * pow(reflected * to_light, mat.getShininess());
							res = res + specular * incoming_light;

						}
					}
				}
			}
			return res * tex_color;
		}

		///////////////////////////////////////////////////////
		// Compute the global illumination of the point
		// - Direct illumination form the emisives surfaces of the scene: OK, but still really noisy
		// - Indirect: Biased
		///////////////////////////////////////////////////////
		RGBColor phong_surface(RayTriangleIntersection const& rti, Ray const& ray, Math::Vector3f & normal, Math::Vector3f & to_view, Math::Vector3f & reflected, size_t depth=0)const
		{
			RGBColor res(0.0, 0.0, 0.0);
			Math::Vector3f const& intersection = rti.intersection();
			Triangle const& tri = *(rti.triangle());
			Material & mat = *(tri.material());

			const double u = rti.uTriangleValue(), v = rti.vTriangleValue();

			RGBColor const tex_color = tri.sampleTexture(u, v);

			const double dist_inter_source = rti.tRayValue();

			//res = res + mat.getAmbient();
			res = res + mat.getEmissive();// / (dist_inter_source + 1);

			//res = res / (dist_inter_source + 1);

			to_view = (ray.source() - intersection).normalized();

			normal = tri.sampleNormal(u, v, ray.source());



			RGBColor mat_diffuse = mat.getDiffuse();
			RGBColor mat_specular = mat.getSpecular();

			reflected = normal * (2 * (normal * to_view)) - to_view;

			if (!(mat_diffuse.isBlack() && mat_specular.isBlack()))
			{
				
				double xi = double(rand()) / double(RAND_MAX);
				if (xi > alpha)
				{
					return res * tex_color;
				}
				//////////////////////////////////////////////
				// Direct illumination (on surface)
				//////////////////////////////////////////////

				for (LightSampler const& m_lightSampler : m_lightSamplers)
				{
					for (unsigned int i = 0; i < m_lightSamples; ++i)
					{
						const Triangle * surface_tri;
						PointLight surface_light;
						m_lightSampler.generate(surface_tri, surface_light);
						Math::Vector3f to_light = (surface_light.position() - rti.intersection());


						if (to_light * normal < 0.f)
						{
							continue;
						}
						double dist_to_light_2 = to_light.norm2(), dist_to_light = sqrt(dist_to_light_2);
						//RGBColor incoming_light = send_ray_to_light(rti, l); // light to tri
						to_light = to_light / dist_to_light;

						RGBColor incoming_light = surface_light.color()* send_ray_to_surface(rti, surface_tri, to_light, dist_to_light_2) / m_lightSamples / alpha;



						if (incoming_light.isBlack())
						{
							continue;
						}

						//Diffuse
						if (!mat_diffuse.isBlack())
						{
							RGBColor diffuse = mat_diffuse * (to_light * normal);
							res = res + diffuse * incoming_light;
						}

						//Specular
						if (!mat_specular.isBlack())
						{

							//maybe compute only the reflected of the ray, instead of the light: only on value to compute, instead on n
							//Math::Vector3f reflected = tri.reflectionDirection(normal, -to_light);
							if (reflected * to_light >= 0)
							{

								///wtf
								RGBColor specular = mat_specular * pow(reflected * to_light, mat.getShininess());
								res = res + specular * incoming_light;

							}
						}
					}
				}

				//////////////////////////////////////////////////////////
				// Indirect illumination
				// TODO 
				// The result is complitly biased
				//////////////////////////////////////////////////////////
				
				
				
				
				Math::RandomDirection diffuse_sampler(reflected);

				Math::Vector3f direction = diffuse_sampler.generate();
				
				
				double diffuse_factor;
				if (direction * reflected < 0)
				{
					direction = -direction;
				}
				diffuse_factor = direction * normal;
				
				if (diffuse_factor < 0)
				{
					diffuse_factor = -diffuse_factor;
					direction = -direction;
				}
				

				res = res + trace_path(Ray(intersection, direction), rti.triangle(), depth+1) * diffuse_factor * mat.getDiffuse() / alpha;
				
			}
			return res * tex_color;
		}





		//////////////////////////////////////////////////////////
		// Phong function, looks into the shadow textures of the triangles
		// The shadows must be pre-compted, else: something like segmentation fault 
		//////////////////////////////////////////////////////////
		RGBColor phong_cache(RayTriangleIntersection const& rti, Ray const& ray, Math::Vector3f & normal, Math::Vector3f & to_view, Math::Vector3f & reflected)const
		{
			RGBColor res(0.0, 0.0, 0.0);
			Math::Vector3f const& intersection = rti.intersection();
			Triangle const& tri = *(rti.triangle());
			Material & mat = *(tri.material());

			double u = rti.uTriangleValue(), v = rti.vTriangleValue();

			RGBColor const tex_color = tri.sampleTexture(u, v);

			double dist_inter_source = rti.tRayValue();

			res = res + mat.getAmbient();
			res = res + mat.getEmissive();

			res = res / (dist_inter_source + 1);

			to_view = (ray.source() - intersection).normalized();

			normal = tri.sampleNormal(u, v, ray.source());



			RGBColor mat_diffuse = mat.getDiffuse();
			RGBColor mat_specular = mat.getSpecular();

			reflected = normal * (2 * (normal * to_view)) - to_view;

			//TriangleMultipleTexture<bool> const& light_texture = tri.m_light_tex;

			if (!(mat_diffuse.isBlack() && mat_specular.isBlack()))
			{

				for (unsigned int i = 0; i < m_lights.size(); ++i)
				{
					PointLight const& l = m_lights[i];
					Math::Vector3f to_light = (l.position() - intersection).normalized();

					if (to_light * normal < 0.f)
					{
						continue;
					}


					double is_lighted = tri.m_light_tex->get_interpolation(tri.get_light_tex_index(i, tri.facing(to_view)), u, v);

					RGBColor incoming_light = l.color() / (l.position() - intersection).norm() * is_lighted;


					if (incoming_light.isBlack())
					{
						continue;
					}

					//Diffuse
					if (!mat_diffuse.isBlack())
					{
						RGBColor diffuse = mat_diffuse * (to_light * normal);
						res = res + diffuse * incoming_light;
					}

					//Specular
					if (!mat_specular.isBlack())
					{

						//maybe compute only the reflected of the ray, instead of the light: only on value to compute, instead on n
						//Math::Vector3f reflected = tri.reflectionDirection(normal, -to_light);
						if (reflected * to_light >= 0)
						{

							///wtf
							RGBColor specular = mat_specular * pow(reflected * to_light, mat.getShininess());
							res = res + specular * incoming_light;

						}
					}
				}
			}
			return res * tex_color;
		}



		///////////////////////////////
		// Computes the number of intersection between the ray and the acceleration structure
		// This version shows the bounding boxes of the geometrys
		///////////////////////////////
		virtual unsigned int count_bb_intersections(Ray const& ray)const
		{
			double t_entry, t_exit;
			unsigned int res = 0;
			for (auto const& p : m_geometries)
			{
				if (p.first.intersect(ray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
				{
					++res;
				}
			}
			return res;
		}


		//////////////////////////////////////
		// counts the bounding boxes of the acceleration structures for the ray, and convets it into a color
		//////////////////////////////////////
		RGBColor count_box(Ray const& ray)const
		{
			int n = 2*count_bb_intersections(ray);

			double red(0), green(0), blue(0);

			if (n < 256)
			{
				blue = n;
			}
			else if (n < 512)
			{
				blue = 512 - n;
				red = n - 256;
			}
			else
			{
				red = 3 * 256 - n;
				green = n - 512;
			}
			

			RGBColor res(red, green, blue);
			res = res / 255. * 2.;
			return res;
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor Scene::sendRay(Ray const & ray, double limit, int depth, int maxDepth)
		///
		/// \brief	Sends a ray in the scene and returns the computed color
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	ray			The ray.
		/// \param	depth   	The current depth.
		/// \param	maxDepth	The maximum depth.
		///
		/// \return	The computed color.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor sendRay(Ray const & ray, int maxDepth, std::stack<const medium*> medium_stack, const Triangle * current_triangle = nullptr, int depth=0)const
		{
			if (depth > maxDepth)
			{
				return RGBColor();
			}
			
			
			
			
			RGBColor result;
			RayTriangleIntersection intersect_info;
			if (intersection(ray, intersect_info, current_triangle))
			{
				Math::Vector3f normal;
				Math::Vector3f to_view;
				Math::Vector3f reflected;

				result = m_phong(this, intersect_info, ray, normal, to_view, reflected);

				Material const& mat = *intersect_info.triangle()->material();
				RGBColor const& spec_color = mat.getSpecular();
				
				if (mat.is_transparent())
				{
					//refraction, but not used at the moment
					medium const& current_medium = *(medium_stack.top());

					medium const& next_medium = *mat.m_medium;

					double this_n = current_medium.index();
					double that_n = next_medium.index();

					std::stack<const medium*> refraction_stack(medium_stack);

					refraction_stack.push(NULL);
					//copy the stack
					Math::Vector3f direction;

					result = result + spec_color * sendRay(Ray(intersect_info.intersection(), ray.invDirection()), maxDepth, medium_stack, intersect_info.triangle(), depth+1);
					

				}
				//*/

				else if (!spec_color.isBlack())
				{
					//increase the stack
					Math::Vector3f reflected = normal * (2 * (normal * to_view)) - to_view;
					Ray reflected_ray(intersect_info.intersection(), reflected);
					result = result + sendRay(reflected_ray, maxDepth, medium_stack, intersect_info.triangle(), depth+1) * spec_color;
				}
			}
			else
			{
				result = get_background_color(ray);
			}
			
			return result;
		}


		double alpha = 0.5;

		////////////////////////////////////////////////////////
		//TODO add refraction
		////////////////////////////////////////////////////////
		RGBColor trace_path(Ray const & ray, const Triangle * current_triangle = nullptr, int depth = 0)const
		{
			//TODO add russian roulette to stop efficiently, and unbiasedly
			if (depth > max_depth)
			{
				return RGBColor();
			}

			


			RGBColor result;
			RayTriangleIntersection intersect_info;
			if (intersection(ray, intersect_info, current_triangle))
			{
				Math::Vector3f normal;
				Math::Vector3f to_view;
				Math::Vector3f reflected;

				result = phong_surface(intersect_info, ray, normal, to_view, reflected, depth);
				if (depth > 1)
				{
					result = result / (1 + intersect_info.tRayValue());
				}
			}
			else
			{
				result = get_background_color(ray);
			}
			
			return result;
		}


		RGBColor sendRayBox(Ray const & ray)const
		{
			return count_box(ray);
		}

		//Used to show the uv of the triangles
		RGBColor sendRay_uv(Ray const& ray)const
		{
			RGBColor result(0, 0, 0);
			RayTriangleIntersection intersect_info;
			if (intersection(ray, intersect_info))
			{
				result[0] = intersect_info.uTriangleValue();
				result[2] = intersect_info.vTriangleValue();
				result[1] = 0.0;
			}
			return result;
		}

		//Used to view the normals of the triangles
		RGBColor sendRay_normal(Ray const& ray)const
		{
			RGBColor result(0, 0, 0);
			RayTriangleIntersection intersect_info;
			if (intersection(ray, intersect_info))
			{
				Math::Vector3f normal = intersect_info.triangle()->sampleNormal(intersect_info.uTriangleValue(), intersect_info.vTriangleValue(), intersect_info.triangle()->facing(-ray.direction()));
				//normal = normal.normalized();
				result[0] = 0.5 + 0.5*normal[0];
				result[1] = 0.5 + 0.5*normal[1];
				result[2] = 0.5 + 0.5*normal[2];
			}
			return result / 2;
		}

		/////////////////////////////////////////
		//Show the "Z buffer"
		// TODO: add normalization
		/////////////////////////////////////////
		RGBColor sendRay_z(Ray const& ray)const
		{
			RGBColor res(0, 0, 0);
			RayTriangleIntersection intersect_info;
			if (intersection(ray, intersect_info))
			{
				double z_value = intersect_info.tRayValue();
				
				

				res[0] = res[1] = res[2] = z_value;
			}
			return res;
		}


	public:
		//get the background color of a ray lost in the skybox
		RGBColor get_background_color(Ray const& ray)const
		{
			if (!m_skybox.isValid() || m_background_color.isBlack())
			{
				return m_background_color;
			}
			Math::Vector3f spc = spherical_coordinate(ray.direction());

			//spc[1] = spc[1] / Math::pi;

			spc[1] = 0.5+-0.5*ray.direction().normalized()[2];

			spc[2] = spc[2] / Math::twoPi;

			Math::Vector2f tex_uv = Math::makeVector(spc[2], spc[1]);

			return m_skybox.pixel(tex_uv)*m_background_color;
		}

		////////////////////////////////////////////////////////////////////////
		//An alternative version of intersection that could be faster for lights
		//That could thearically be a little bit faster
		////////////////////////////////////////////////////////////////////////
		virtual bool intersection_light(Ray const& ray, double light_t, const Triangle * current = nullptr)const
		{
			CastedRay cray(ray);
			double current_min_t = std::numeric_limits<double>::max();
			for (const std::pair< BoundingBox, Geometry > & geom : m_geometries)
			{
				double t1, t2;
				if (geom.first.intersect(cray, 0, current_min_t, t1, t2))
				{
					//current_min_t = ::std::min(t1, std::min(t2, current_min_t));

					//geom.second.intersection(cray);
					for (Triangle const& t : geom.second.getTriangles())
					{
						if (&t != current)
						{
							cray.intersect(&t);
						}
					}
					if (cray.validIntersectionFound())
					{
						RayTriangleIntersection const& intersect_info = cray.intersectionFound();
						if (intersect_info.tRayValue() < light_t)
							return true;

						//sale
						current_min_t = std::max(t1, t2);
					}

				}



			}
			return false;
		}



		//////////////////////////////////////
		//Sends a ray to a points light and returns the contribution of the light
		//depending of the visibility of the light
		//////////////////////////////////////
		RGBColor send_ray_to_light(Ray const & ray, PointLight const& l, const Triangle * current_triangle=nullptr)const
		{
			//return RGBColor(1.0, 1.0, 1.0);
			double dist_to_light_2 = (l.position() - ray.source()).norm2();
			double dist_to_light = sqrt(dist_to_light_2);

			RGBColor res = l.color() / dist_to_light;


			//if (res.grey() < 0.01)//The contribution of light is too small to be calculated...
			//	return RGBColor();

			if (intersection_light(ray, dist_to_light, current_triangle))
			{
				return RGBColor();
			}
			else
			{
				// light is reached
				return res;
			}
		}

		/////////////////////////////////////////////////////////
		//This one sends a ray from the light to the intersection
		// DEPRECATED
		////////////////////////////////////////////////////////
		RGBColor send_ray_from_light(RayTriangleIntersection const& current_inter, PointLight const& l)const
		{
			
			Ray ray(l.position(), current_inter.intersection() - l.position());
			
			RayTriangleIntersection intersect_info;
			intersection(ray, intersect_info);
			
			if (current_inter.triangle() == intersect_info.triangle())
			{
				double dist_to_tight_2 = (l.position() - current_inter.intersection()).norm2();
				return l.color() / sqrt(dist_to_tight_2);
			}
			else
			{
				return RGBColor();
			}
		}


		////////////////////////////////////////////
		//Sends a ray to a sampled point of an emissive surface
		////////////////////////////////////////////
		RGBColor send_ray_to_surface(RayTriangleIntersection const& current_inter, const Triangle * surface_tri, const Math::Vector3f & to_light, double dist_to_light_2)const
		{
			
			Ray ray(current_inter.intersection(), to_light);
			RayTriangleIntersection intersect_info;
			intersection(ray, intersect_info);
			if (intersect_info.valid() && surface_tri == intersect_info.triangle())
			{
				return surface_tri->material()->getEmissive() / intersect_info.tRayValue();
			}
			return 0;
		}




		size_t num_shadow = 0;



	protected:


		//////////////////////////////////////////////
		//Computes the visibily between a point and a light
		//////////////////////////////////////////////
		bool is_lighted(Math::Vector3f const& point, PointLight const& l, Triangle * current_triangle = nullptr)const
		{
			Math::Vector3f to_light = l.position() - point;

			double dist_to_light_2 = to_light.norm2();
			double dist_to_light = sqrt(dist_to_light_2);


			Ray ray(point, to_light / dist_to_light);

			

			if (intersection_light(ray, dist_to_light, current_triangle))
			{
				return 0;
			}
			else
			{
				// light is reached
				return 1;
			}
		}


		

		/////////////////////////////////////////////////
		//Compute the shadow texture of a triangle,
		//TODO:
		// - Optimisation (allocate if necessary, use more cache efficient representation)
		// - CLean the artefacts
		////////////////////////////////////////////////
		void compute_triangle_shadow(Triangle & triangle, size_t width, size_t height)
		{
			if (triangle.m_light_tex != nullptr)
			{
				delete triangle.m_light_tex;
			}
			triangle.m_light_tex = new TriangleMultipleTexture<bool>(2 * m_lights.size(), width, height);
			TriangleMultipleTexture<bool> & light_texture = *(triangle.m_light_tex);

			for (unsigned char facing = 0; facing < 2; ++facing)
			{
				
				for (unsigned int light_index = 0; light_index < m_lights.size(); ++light_index)
				{
					PointLight const& l = m_lights[light_index];
					bool facing_light = false;
					for (unsigned int i = 0; i < 3; ++i)
					{
						Math::Vector3f vertex = triangle.vertex(i);
						Math::Vector3f normal = triangle.getVertexNormal(i, facing);

						Math::Vector3f to_light = l.position() - vertex;

						if (to_light * normal > 0.0)
						{
							facing_light = true;
							break;
						}
					}

					size_t tex_index = triangle.get_light_tex_index(light_index, facing);

					
					const double mid_pixel = 0.;
					
					
					if (facing_light)
					{
						for (size_t v = 0; v < light_texture.height(); ++v)
						{
							double v_value = light_texture.get_v(v + mid_pixel);
							for (size_t u = 0; u < light_texture.width(tex_index, v); ++u)
							{
								
								
								double u_value = light_texture.get_u(u + mid_pixel);
								Math::Vector3f sampled_point = triangle.samplePoint(u_value, v_value);
								Math::Vector3f sampled_normal = triangle.sampleNormal(u_value, v_value, facing);

								Math::Vector3f to_light = l.position() - sampled_point;

								if (to_light * sampled_normal < 0.0)
								{
									light_texture.set_pixel(tex_index, u, v, false);
								}
								else
								{

									light_texture.set_pixel(tex_index, u, v, is_lighted(sampled_point, l, &triangle));
								}
							}
						}
						if (light_texture.set_uniform(tex_index))
						{
							
						}
						else
						{
							num_shadow += 1;

						}
					}
					else
					{
						light_texture.clear(tex_index);
						light_texture.create_default_values();
						light_texture.set_default(tex_index, false);
					}
					
				}
			}
		}

	public:



		/////////////////////////////////////
		//Pre compute the shadows of all the triangles of the scene
		// - This function can be quite expenseive
		// - The size of the shades are at the scales width_factor and height_factor
		// TODO: add a minimal size
		////////////////////////////////////
		void pre_compute_shadows(precision width_factor, precision height_factor, size_t min_width=1, size_t min_height=1)
		{
			num_shadow = 0;
			std::vector<Triangle *> triangle_list;
			triangle_list.reserve(m_number_triangle);
			for (auto & p : m_geometries)
			{
				for (Triangle & triangle : p.second.getTriangles())
				{
					triangle_list.push_back(&triangle);
				}
			}

			int step = std::max(1ull, triangle_list.size() / 100ull);

#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < triangle_list.size(); ++i)
			{
				Triangle & triangle = *(triangle_list[i]);

				size_t width = std::max(min_width, (size_t)(width_factor * triangle.uAxis().norm()) + 1);
				size_t height = std::max(min_height, (size_t)(height_factor * triangle.vAxis().norm()) + 1);

				compute_triangle_shadow(triangle, width, height);

				//Printing status
				if (i % step == 0)
				{
					int percent = (i * 100) / triangle_list.size();
#pragma omp critical(cout)
					std::cout << "\r" + progession_bar(i, triangle_list.size(), 50) << std::flush;
				}
			}
			std::cout << "\r" + progession_bar(1, 1, 50) << std::endl;
			m_phong = &Scene::phong_cache;
		}




#define PHONG_CACHE 0
#define PHONG_RAY 1
#define PHONG_SURFACE 2


		//TODO use an array
		//Switches the phong fonction to be used, will soon be done better
		void set_phong(unsigned char mode)
		{
			switch (mode)
			{
			case PHONG_CACHE:
				m_phong = &Scene::phong_cache;
				break;
			case PHONG_RAY:
				m_phong = &Scene::phong_ray;
				break;
			case PHONG_SURFACE:
				//m_phong = &Scene::phong_surface;
			default:
				break;
			}
		}


		


		

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::compute(int maxDepth)
		///
		/// \brief	Computes a rendering of the current scene, viewed by the camera.
		/// 		
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	maxDepth	The maximum recursive depth.
		/// \param  subPixelDivision subpixel subdivisions to handle antialiasing
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void compute(int maxDepth, unsigned char render_mode=1, int subPixelDivision = 1, int passPerPixel = 1)
		{
			// We prepare the light sampler (the sampler only stores triangles with a non null emissive component).
			for (auto it = m_geometries.begin(), end = m_geometries.end(); it != end; ++it)
			{
				LightSampler ls;
				ls.add(it->second);
				if (ls.hasLights())
				{
					m_lightSamplers.push_back(ls);
				}
			}

			//set_phong(PHONG_SURFACE);

			// Step on x and y for subpixel sampling
			double step = 1.0f/subPixelDivision ;
			// Table accumulating values computed per pixel (enable rendering of each pass)
			::std::vector<::std::vector<::std::pair<int, RGBColor> > > pixelTable(m_visu->width(), ::std::vector<::std::pair<int, RGBColor> >(m_visu->width(), ::std::make_pair(0, RGBColor()))) ;

			// 1 - Rendering time
			LARGE_INTEGER frequency;        // ticks per second
			LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			QueryPerformanceFrequency(&frequency);
			// start timer
			QueryPerformanceCounter(&t1);
			// Rendering pass number
			m_pass = 0;
			// Rendering

			std::stack<const medium*> medium_stack;
			medium_stack.push(&m_medium);

			bool enter = false;

			for(int passPerPixelCounter = 0 ; passPerPixelCounter<passPerPixel ; ++passPerPixelCounter)
			{
				for (double xp = -0.5; xp < 0.5; xp += step)
				{
					for (double yp = -0.5; yp < 0.5; yp += step)
					{
						::std::cout << "Pass: " << m_pass << "/" << passPerPixel * subPixelDivision * subPixelDivision << ::std::endl;
						++m_pass;
						// Sends primary rays for each pixel (uncomment the pragma to parallelize rendering)
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
						for (int y = 0; y < m_visu->height(); y++)
						{
							for (int x = 0; x < m_visu->width(); x++)
							{
//#pragma omp critical (visu)
								//m_visu->plot(x*scale, y*scale, RGBColor(1000.0, 0.0, 0.0), scale);
								// Ray casting
								RGBColor result;
								if (render_mode == 1)
								{
									//result = sendRay(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()), maxDepth, medium_stack);
									result = trace_path(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()));
								}
								else if(render_mode == 0)
								{
									result = sendRayBox(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()));
								}
								else if (render_mode == 2)
								{
									result = sendRay_uv(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()));
								}
								else if (render_mode == 3)
								{
									result = sendRay_normal(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()));
								}
								else if (render_mode == 4)
								{
									result = sendRay_z(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()));
								}
								
								//std::cout << "Result is black: " << result.isBlack() << std::endl;
								// Accumulation of ray casting result in the associated pixel
								
								//result = result * result;
								
								::std::pair<int, RGBColor> & currentPixel = pixelTable[x][y];
								currentPixel.first++;
								currentPixel.second = currentPixel.second + result;
								// Pixel rendering (with simple tone mapping)
#pragma omp critical (visu)
								m_visu->plot(x, y, pixelTable[x][y].second / (double)(pixelTable[x][y].first) * 10);
								// Updates the rendering context (per pixel) - warning per pixel update can be costly...
//#pragma omp critical (visu)
								//m_visu->update();
							}
							// Updates the rendering context (per line)
//#pragma omp critical (visu)
//							m_visu->update();
						}
						// Updates the rendering context (per pass)

						enter = m_visu->update();
						// We print time for each pass
						QueryPerformanceCounter(&t2);
						elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
						double remainingTime = (elapsedTime / m_pass)*(passPerPixel * subPixelDivision * subPixelDivision - m_pass);
						::std::cout << "time: " << elapsedTime << "s. " <<", remaining time: "<< remainingTime << "s. " <<", total time: "<< elapsedTime + remainingTime << ::std::endl;

						if (enter)
						{
							goto compute__end_loop;
						}
						
					}
					
				}
				
			}

		compute__end_loop:


			// stop timer
			QueryPerformanceCounter(&t2);
			elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
			::std::cout<<"time: "<<elapsedTime<<"s. "<<::std::endl ;

			while (!enter)
			{
				enter = m_visu->update();
			}

			//set_phong(PHONG_RAY);

			return;
		}







		void realtime_compute(int maxDepth = 1, unsigned char render_mode=1)
		{
			if (render_mode == 0) // box mode
			{
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
				for (int y = 0; y < m_visu->height(); y++)
				{
					for (int x = 0; x < m_visu->width(); x++)
					{
						RGBColor result = sendRayBox(m_camera.getRay(((double)x) / m_visu->width(), ((double)y) / m_visu->height()));


#pragma omp critical (visu)
						m_visu->plot(x, y, result * 10);

					}

				}
			}
			else if (render_mode == 1) // classic mode
			{
				std::stack<const medium*> medium_stack;
				medium_stack.push(&m_medium);
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
				for (int y = 0; y < m_visu->height(); y++)
				{
					for (int x = 0; x < m_visu->width(); x++)
					{
						RGBColor result = sendRay(m_camera.getRay(((double)x) / m_visu->width(), ((double)y) / m_visu->height()), maxDepth, medium_stack);


#pragma omp critical (visu)
						m_visu->plot(x, y, result * 10);

					}

				}
			}
			else if (render_mode == 2) // uv mode
			{
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
				for (int y = 0; y < m_visu->height(); y++)
				{
					for (int x = 0; x < m_visu->width(); x++)
					{
						RGBColor result = sendRay_uv(m_camera.getRay(((double)x) / m_visu->width(), ((double)y) / m_visu->height()));


#pragma omp critical (visu)
						m_visu->plot(x, y, result * 10);

					}

				}
			}
			else if (render_mode == 3) // normal mode
			{
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
				for (int y = 0; y < m_visu->height(); y++)
				{
					for (int x = 0; x < m_visu->width(); x++)
					{
						RGBColor result = sendRay_normal(m_camera.getRay(((double)x) / m_visu->width(), ((double)y) / m_visu->height()));


#pragma omp critical (visu)
						m_visu->plot(x, y, result * 10);

					}

				}
			}
			else if (render_mode == 4) // z buffer
			{
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
				for (int y = 0; y < m_visu->height(); y++)
				{
					for (int x = 0; x < m_visu->width(); x++)
					{
						RGBColor result = sendRay_z(m_camera.getRay(((double)x) / m_visu->width(), ((double)y) / m_visu->height()));


#pragma omp critical (visu)
						m_visu->plot(x, y, result * 10);

					}

				}
			}

			m_visu->update();
				

			m_lightSamplers.clear();
			
		}

	} ;
}

#endif
