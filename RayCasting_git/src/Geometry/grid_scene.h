#pragma once

#include "Scene.h"
#include <utils.h>
#include <set>
#include <limits>



namespace Geometry
{
	class grid_scene : public Scene
	{
	public:

		grid_scene(Visualizer::Visualizer * visu):
			Scene(visu), grid(nullptr)
		{}

		virtual ~grid_scene()
		{
			assert(grid != nullptr);
			for (unsigned int i = 0; i < grid_res[0]; ++i)
			{
				for (unsigned int j = 0; j < grid_res[1]; ++j)
					delete[] grid[i][j];
				delete[] grid[i];
			}
			delete[] grid;
		}

		


		virtual bool intersection(Ray const& ray, RayTriangleIntersection & intersect_info, Triangle * current = nullptr)
		{
			//check if camera ouside the box
			CastedRay cray(ray);

			Math::Vector3f dv = cray.direction().normalized();

			Math::Vector3f current_pos = cray.source() - min;

			Math::Vector<unsigned int, 3> current_cell = get_current_cell(current_pos);

			Math::Vector3f ts;

			Math::Vector3f dt;

			double current_t = std::numeric_limits<double>::max();


			for (int i = 0; i < 3; ++i)
			{
				
				if (dv[i] < 0.f)
				{
					ts[i] = (current_pos[i] - current_cell[i] * cell_size[i]) / (-dv[i]);
					dt[i] = cell_size[i] / (-dv[i]);
				}
				else if (dv[i] == 0.0)
				{
					ts[i] = std::numeric_limits<double>::max();
					dt[i] = std::numeric_limits<double>::max();
				}
				else
				{
					dt[i] = cell_size[i] / dv[i];
					ts[i] = (cell_size[i] - (current_pos[i] - current_cell[i] * cell_size[i])) / (dv[i]);
				}
			}

			current_t = 0;
			bool valid_intersect_found = false;
			double valid_t = std::numeric_limits<double>::max();

			while (true)
			{
				current_t = std::min(ts[0], std::min(ts[1], ts[2]));

				if (cray.validIntersectionFound())
				{
					intersect_info = cray.intersectionFound();
					valid_t = (ray.source() - intersect_info.intersection()).norm();
					if (valid_t > current_t)
					{
						//std::cout << true << std::endl;
						return true;
					}
				}

				if (current_cell[0] < 0 || current_cell[0] >= grid_res[0] || current_cell[1] < 0 || current_cell[1] >= grid_res[1] || current_cell[2] < 0 || current_cell[2] >= grid_res[2])
				{
					//out of the box
					return false;
				}

				std::vector<const Triangle*> & cell_content = grid[current_cell[0]][current_cell[1]][current_cell[2]];
				for (const Triangle * tri : cell_content)
				{
					cray.intersect(tri);
				}
				
				
				//move to the next cell
				unsigned int min_index = 0;
				double min_t = ts[min_index];
				for (unsigned int i = 1; i < 3; ++i)
				{
					if (ts[i] < min_t)
					{
						min_index = 0;
						min_t = ts[min_index];
					}
				}

				if (dv[min_index] < 0)
				{
					--current_cell[min_index];
				}
				else
				{
					++current_cell[min_index];
				}
					
				//update ts[i]
				ts[min_index] += dt[min_index];
				
			}
			
		}


		void pre_compute(float alpha=4)
		{
			/////////////////////////////////////////////////////////////
			//Pardon
			/*
			for (auto & p : m_geometries)
			{
				p.second.translate(-m_sceneBoundingBox.min());
				p.first = BoundingBox(p.second);
			}
			std::vector<PointLight> old_lights = m_lights;
			m_lights.clear();
			for (PointLight & l : old_lights)
			{
				m_lights.push_back(PointLight(l.position() - m_sceneBoundingBox.min(), l.color()));
			}
			m_camera.setPosition(m_camera.getPosition() - m_sceneBoundingBox.min());

			m_sceneBoundingBox = BoundingBox(Math::makeVector(0.f, 0.f, 0.f), m_sceneBoundingBox.max() - m_sceneBoundingBox.min());
			*/
			
			/////////////////////////////////////////////////////////////

			assert(grid == nullptr);

			min = m_sceneBoundingBox.min();

			unsigned int triangle_number = 0;
			unsigned int total_pointers = 0;
			for (auto & p : m_geometries)
			{
				triangle_number += p.second.getTriangles().size();
			}
			Math::Vector3f box_size = m_sceneBoundingBox.max() - min;
			std::cout << "Box size: " << box_size << std::endl;
			double box_volume = box_size[0] * box_size[1] * box_size[2];

			Math::Vector3f grid_float_res = box_size * sqrt(alpha * triangle_number / box_volume);
			grid_res = Math::makeVector<unsigned int>((unsigned int)(grid_float_res[0] + .5), (unsigned int)(grid_float_res[1] + .5), (unsigned int)(grid_float_res[2] + .5));

			std::cout << "Resolution of the grid: " << grid_res << std::endl;

			cell_size = Math::makeVector<double>(box_size[0] / grid_res[0], box_size[1] / grid_res[1], box_size[2] / grid_res[2]);
			std::cout << "Nuber of cells: " << grid_res[0] * grid_res[1] * grid_res[2] << std::endl;
			// Maybe have contiguous data
			//tri dimentional array allocation 
			grid = new std::vector<const Triangle*>**[grid_res[0]];
			for (unsigned int i = 0; i < grid_res[0]; ++i)
			{
				grid[i] = new std::vector<const Triangle*>*[grid_res[1]];
				for (unsigned int j = 0; j < grid_res[1]; ++j)
				{
					grid[i][j] = new std::vector<const Triangle*>[grid_res[2]];
				}
			}
			
			

			for (auto & p : m_geometries)
			{
				std::cout << "Geometry: " << &p.second <<", ("<< p.second.getTriangles().size()<< " triangles)" << std::endl;
				for (const Triangle & tri : p.second.getTriangles())
				{
					//std::cout << "Triangle " << &tri << ": ";
					//std::cout << "(" << tri.vertex(0) << ", " << tri.vertex(1) << ", " << tri.vertex(2) << ")"<< "center: "<<tri.center()<<std::endl;
					for (unsigned int i = 0; i < grid_res[0]; ++i)
					{
						for (unsigned int j = 0; j < grid_res[1]; ++j)
						{
							for (unsigned int k = 0; k < grid_res[2]; ++k)
							{
								Math::Vector3f base = min + Math::makeVector(i*cell_size[0], j*cell_size[1], k*cell_size[2]);
								//std::cout << base << std::endl;
								BoundingBox bb(base, base + cell_size);
								if (collide(bb, tri))
								{
									
									//std::cout<< " to the cell (" << i << ", " << j << ", " << k << ")" << std::endl;
									++total_pointers;
									grid[i][j][k].push_back(&tri);
								}
							}
						}
					}
				}
			}


			/*
			for (auto& p : m_geometries)
			{
				for (Triangle const& t : p.second.getTriangles())
				{
					std::cout << "Triangle: " << &t << std::endl;
					BoundingBox bb = get_bounding_box(t);
					std::vector<Math::Vector3f> vertices;
					Math::Vector3f tab[2] = { bb.min() - min, bb.max() - min};//Les 8 sommets de la bb
					for (unsigned char i = 0; i < 2; ++i)
						for (unsigned char j = 0; j < 2; ++j)
							for (unsigned char k = 0; k < 2; ++k)
								vertices.push_back(Math::makeVector<float>(tab[i][0], tab[j][1], tab[k][2]));

					std::set<std::vector<const Triangle*>*> containing_grid;
					// faut faire un collision de bb ... ou pas (boucles sur les inters ?? dans la grid) 
					for (Math::Vector3f & vertex : vertices)
					{
						Math::Vector<unsigned int, 3> vertex_grid_pos = get_current_cell(vertex);
						
						// Turbo pas sur
						std::cout << "vertex_grid_pos: " << vertex_grid_pos << std::endl;
						std::vector<const Triangle*>* pointer = &(grid[vertex_grid_pos[0]][vertex_grid_pos[1]][vertex_grid_pos[2]]);
						containing_grid.insert(pointer);
						
					}
					total_pointers += containing_grid.size();
					for (std::vector<const Triangle*> * vec : containing_grid)
					{
						vec->push_back(&t);
					}
				}
			}
			*/
			
			std::cout << "Pre compute Done!" << std::endl;
			std::cout << "Mean number of pointer per triangle: " << ((float)total_pointers) / ((float)triangle_number) << std::endl;

		}

	protected:

		Math::Vector3f min;
		Math::Vector<unsigned int, 3> grid_res;
		Math::Vector3f cell_size;
		std::vector<const Triangle*> *** grid;


		Math::Vector<unsigned int, 3> get_current_cell(Math::Vector3f vertex)
		{
			Math::Vector<unsigned int, 3> res;
			for (int i = 0; i < 3; ++i)
			{
				res[i] = (unsigned int)(vertex[i] / cell_size[i]);
				//sale
				if (res[i] == grid_res[i])
				{
					--res[i];
				}
			}
			return res;
		}



		bool collide(BoundingBox const& bb, Triangle const& t)
		{
			double bb_radius2 = ((bb.max() - bb.min())/2).norm2();
			double tri_rad2 = 0;
			for (unsigned int i=0; i<3; ++i)
			{
				Math::Vector3f const& vertex = t.vertex(i);
				double dist = (vertex - t.center()).norm2();
				if (dist > tri_rad2)
				{
					tri_rad2 = dist;
				}
			}

			double dist_between_objects = (((bb.max() + bb.min()) / 2) - t.center()).norm2();
			if (dist_between_objects > bb_radius2 + tri_rad2)
			{
				return false;
			}


			return true;
		}
	
	};
}