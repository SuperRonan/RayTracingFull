#pragma once


#include "Scene.h"
#include <vector>
#include "../struct/static_dual_tree.h"
#include <functional>
#include <deque>
#include <algorithm>
#include <climits>
#include <string>

namespace Geometry 
{

	typedef std::vector<const Triangle*> primitive_collection;

	/*
	typedef static_dual_tree<2, BoundingBox, primitive_collection> bvh_type;
	typedef static_dual_node<2, BoundingBox, primitive_collection> bvh_node;
	typedef static_dual_leaf<2, BoundingBox, primitive_collection> bvh_leaf;
	*/
	typedef turbo_static_dual_tree<2, BoundingBox, primitive_collection> bvh_type;



	///////////////////////////////////////////////////////
	// A scene with a BVH acceleration structure
	//////////////////////////////////////////////////////
	class bvh_scene : public Scene
	{
	public:

		

		bvh_scene(Visualizer::Visualizer * visu) :
			Scene(visu), bvh(nullptr)
		{}

		virtual ~bvh_scene()
		{
			if (bvh != nullptr)
			{
				delete bvh;
			}
		}


		void printStats()
		{
			//maybe print the (basic) struct of the tree...
			::std::cout << "Scene bvh: " << m_number_triangle << " triangles" << ::std::endl;
		}



		//
		virtual bool intersection(Ray const& ray, RayTriangleIntersection & intersect_info, const Triangle * current = nullptr)const
		{
			CastedRay cray(ray);
			double t_entry;
			double t_exit;
			
			if (m_sceneBoundingBox.intersect(cray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
			{
				intersection(cray, current, bvh, t_entry, t_exit);
				if (cray.validIntersectionFound())
				{
					intersect_info = cray.intersectionFound();
					return true;
				}
			}
			return false;
			
		}

	protected:

		///////////////////////////////////////////////////
		// Compute the intersection of the ray with the scene par going through the BVH
		//////////////////////////////////////////////////
		void intersection(CastedRay & cray, const Triangle * current_tri, const bvh_type const * current_node, double t_entry, double t_exit)const
		{
			//std::cout << "depth: " << depth << std::endl;
			if (current_node->is_leaf())
			{
				//auto s = std::to_string(depth) + "|";
				//std::cout << s<<std::endl;
				primitive_collection const& triangles = current_node->get_leaf_value();
				for (const Triangle* tri : triangles)
				{
					//if (tri->normal() * cray.direction() < 0.0)
						if(tri != current_tri)
							cray.intersect(tri);
				}
			}
			else
			{
				double t_entry_l, t_exit_l, t_entry_r, t_exit_r;
				//compute ray intersection with the sons of the node

				//bb intersections

				bool intersect_l = (*current_node)[0]->get_node_value().intersect(cray, t_entry, t_exit, t_entry_l, t_exit_l);
				bool intersect_r = (*current_node)[1]->get_node_value().intersect(cray, t_entry, t_exit, t_entry_r, t_exit_r);

				if (intersect_l)
				{
					if (intersect_r)//both
					{
						//maybe go to right instead of left sometimes... dependings on the t_entrys

						//begin by the "closest box"
						if (t_entry_l < t_entry_r) //begin by the left box 
						{
							intersection(cray, current_tri, (*current_node)[0], t_entry_l, t_exit_l);

							if (cray.validIntersectionFound())
							{
								RayTriangleIntersection const& inter = cray.intersectionFound(); //aka inter_l
								if (inter.tRayValue() > t_entry_r)
								{
									intersection(cray, current_tri, (*current_node)[1], t_entry_r, t_exit_r);
								}
							}
							else // no real left intersection 
							{
								intersection(cray, current_tri, (*current_node)[1], t_entry_r, t_exit_r);
							}
						}
						else // begin by the right box
						{
							intersection(cray, current_tri, (*current_node)[1], t_entry_r, t_exit_r);
							if (cray.validIntersectionFound())
							{
								const RayTriangleIntersection & inter = cray.intersectionFound();
								if (inter.tRayValue() > t_entry_l)
								{
									intersection(cray, current_tri, (*current_node)[0], t_entry_l, t_exit_l);
								}
							}
							else
							{
								intersection(cray, current_tri, (*current_node)[0], t_entry_l, t_exit_l);
							}
						}
					}
					else//only l
					{
						intersection(cray, current_tri, (*current_node)[0], t_entry_l, t_exit_l);
					}
				}
				else
				{
					if (intersect_r)//only r
					{
						intersection(cray, current_tri, (*current_node)[1], t_entry_r, t_exit_r);
					}
					else//no intersection
					{
						return;
					}
				}
			}
		}


		virtual bool intersection_light(Ray const& ray, double light_t, const Triangle * current = nullptr)const
		{
			CastedRay cray(ray);
			double t_entry;
			double t_exit;

			if (m_sceneBoundingBox.intersect(cray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
			{
				t_exit = std::min(t_exit, light_t);
				intersection(cray, current, bvh, t_entry, t_exit);
				if (cray.validIntersectionFound())
				{
					//assert(light_t <= cray.intersectionFound().tRayValue());
					if(cray.intersectionFound().tRayValue() < light_t)
						return true;
				}
			}
			return false;
		}


		bool intersect(double a, double b, double c, double d)
		{
			return a <= d && b >= c;
		}


	public:

		unsigned int count_bb(Ray const& ray, const bvh_type const * node)const
		{
			unsigned int res = 0;
			double t0, t1;
			if (node->get_node_value().intersect(ray, 0, std::numeric_limits<double>::max(), t0, t1))
			{
				res = 1;
				if (node->is_node())
				{
					res += count_bb(ray, (*node)[0]);
					res += count_bb(ray, (*node)[1]);
				}
				
			}
			return res;
		}

		virtual unsigned int count_bb_intersections(const Ray & ray)const
		{
			return count_bb(ray, bvh);
		}

		
		////////////////////////////////
		// Computes the BVH
		// The parameters are for the SAH
		////////////////////////////////
		void pre_compute(int n_max, int ct, int ci)
		{
			primitive_collection current_primitives;
			for (auto & p : m_geometries)
			{
				//current_primitives.insert(current_primitives.end(), p.second.getTriangles().begin(), p.second.getTriangles().end());
				for (const Triangle & tri : p.second.getTriangles())
				{
					current_primitives.push_back(&tri);
				}
			}

			assert(bvh == nullptr);
			assert(current_primitives.size() >= 1);


			pre_compute_body(n_max, ct, ci, current_primitives, m_sceneBoundingBox, bvh, 0);

			
			

		}

	protected:

		
		double best_sah_partition(int ci, int ct, primitive_collection const& sorted_primitives, double current_box_surface, primitive_collection & p_left, primitive_collection & p_right, BoundingBox & bb_left, BoundingBox & bb_right)
		{
			assert(sorted_primitives.size() > 1);
			//assert(p_left.empty());
			//assert(p_right.empty());

			std::deque<const Triangle *> tmp_left, tmp_right;

			

			double best_sah = std::numeric_limits<double>::max();

			BoundingBox * right_bb = new BoundingBox[sorted_primitives.size()];

			BoundingBox tmp;
			
			for (int i = sorted_primitives.size() - 1; i >= 0; --i)
			{
				tmp_right.push_front(sorted_primitives[i]);
				tmp.update(*(sorted_primitives[i]));
				right_bb[i] = tmp;
			}

			

			BoundingBox left_bb;

			BoundingBox best_left, best_right;

			unsigned int best_k(0);

			for (unsigned int k = 1; k < sorted_primitives.size(); ++k)
			{
				tmp_left.push_back(tmp_right.front());
				left_bb.update(*tmp_right.front());
				tmp_right.pop_front();

				
				double left_bb_surface = left_bb.surface();
				double right_bb_surface = right_bb[k].surface();

				double current_sah = sah(ct, ci, tmp_left.size(), tmp_right.size(), current_box_surface, left_bb_surface, right_bb_surface);

				if(current_sah < best_sah)
				{
					best_sah = current_sah;
					best_left = left_bb;
					best_k = k;
				}
				else
				{
					
				}
			}
			tmp_left.clear();
			tmp_right.clear();
			
			best_right = right_bb[best_k];
			for (unsigned int i = 0; i < best_k; ++i)
			{
				p_left.push_back(sorted_primitives[i]);
				//best_left.update(*sorted_primitives[i]);
			}
			for (unsigned int i = best_k; i < sorted_primitives.size(); ++i)
			{
				p_right.push_back(sorted_primitives[i]);
			}
	
			bb_left = best_left;
			bb_right = best_right;

			//assert
			if (bb_left != construct_bb(p_left))
			{
				std::cout << "left not OK" << std::endl;
			}

			if (bb_right != construct_bb(p_right))
			{
				std::cout << "right not OK" << std::endl;
			}

			delete[] right_bb;

			return best_sah;
		}


		double sah(int ct, int ci, unsigned int nl, unsigned int nr, double sb, double sbl, double sbr)const
		{
			return ct + ci * (sbl / sb * nl + sbr / sb * nr);
		}


		/*
		std::function<bool (const Triangle*, const Triangle *)> comp[] = { (const Triangle* a, const Triangle * b) {
			return a->center()[0] < b->center()[0];
		} };
		*/

		typedef bool(*turbo_comparator) (const Triangle*, const Triangle*);

		static bool comp_0(const Triangle* a, const Triangle* b)
		{
			return a->center()[0] < b->center()[0];
		}

		static bool comp_1(const Triangle* a, const Triangle* b)
		{
			return a->center()[1] < b->center()[1];
		}

		static bool comp_2(const Triangle* a, const Triangle* b)
		{
			return a->center()[2] < b->center()[2];
		}

		turbo_comparator comparators[3] = {comp_0, comp_1, comp_2};



		BoundingBox construct_bb(primitive_collection const& pc)
		{
			BoundingBox res;
			for (const Triangle* t : pc)
			{
				res.update(*t);
			}
			return res;
		}

		

		/////////////////////////////////
		//Finds the best partition of collection of primitives using the SAH
		/////////////////////////////////
		void find_best_partition_with_sah(int ci, int ct, primitive_collection const& current_primitives, 
			BoundingBox const& current_bb, primitive_collection & p_left, primitive_collection & p_right, BoundingBox & bb_left, BoundingBox & bb_right)
		{

			double best_sah = std::numeric_limits<double>::max();

			double current_bb_surface = current_bb.surface();


			for (unsigned int axis = 0; axis < 3; ++axis)//for each axis
			{
				//sort the primitve along the current axis

				primitive_collection axis_sorted = current_primitives;

				std::sort(axis_sorted.begin(), axis_sorted.end(), comparators[axis]);


				//compute the min sah for this axis

				primitive_collection current_p_left, current_p_right;
				BoundingBox current_bb_left, current_bb_right;

				
				double current_sah = best_sah_partition(ci, ct, axis_sorted, current_bb_surface, current_p_left, current_p_right, current_bb_left, current_bb_right);
				
				
				if (current_sah < best_sah)
				{
					best_sah = current_sah;
					p_left = current_p_left;
					p_right = current_p_right;
					bb_left = current_bb_left;
					bb_right = current_bb_right;
				}

			}
		}

		///////////////////////////////////////////
		//Medium partition, sorted aling a random axis
		//Faster to build than the SAH, but less efficient
		///////////////////////////////////////////
		void find_med_partition(primitive_collection const& current_primitives, primitive_collection & p_left, primitive_collection & p_right, 
			BoundingBox & bb_left, BoundingBox & bb_right)
		{
			unsigned int axis = rand() % 3;//L O L
			primitive_collection axis_sorted = current_primitives;

			std::sort(axis_sorted.begin(), axis_sorted.end(), comparators[axis]);

			for (unsigned int i = 0; i < axis_sorted.size() / 2; ++i)//Left
			{
				p_left.push_back(axis_sorted[i]);
				bb_left.update(*axis_sorted[i]);
			}
			for (unsigned int i = axis_sorted.size() / 2; i < axis_sorted.size(); ++i)//Right
			{
				p_right.push_back(axis_sorted[i]);
				bb_right.update(*axis_sorted[i]);
			}
			
		}

		
		virtual void pre_compute_body(int n_max, int ct, int ci, primitive_collection & current_primitives,
			BoundingBox const& current_bb, bvh_type *& current_tree_pos, unsigned int depth)
		{
			assert(current_tree_pos == nullptr);
			assert(current_primitives.size() > 1);
			if (current_primitives.size() <= n_max)
			{

				current_tree_pos = bvh_type::make_new_leaf(current_bb, current_primitives);
			}
			else
			{
				bvh_type * current_node= bvh_type::make_new_node(current_bb);
				current_tree_pos = current_node;

				//sparate the current_primitives:

				primitive_collection p_left, p_right;
				BoundingBox bb_left, bb_right;
				
				find_best_partition_with_sah(ci, ct, current_primitives, current_bb, p_left, p_right, bb_left, bb_right);
				//find_med_partition(current_primitives, p_left, p_right, bb_left, bb_right);

				current_primitives.clear();
				pre_compute_body(n_max, ct, ci, p_left, bb_left, (*current_node)[0], depth+1);
				pre_compute_body(n_max, ct, ci, p_right, bb_right, (*current_node)[1], depth+1);

			}
		}



		bvh_type * bvh;

		
	};
}