#pragma once



#include "bvh_scene.h"

#include "tbb/parallel_invoke.h"


namespace Geometry
{

	//////////////////////////////////////////////////
	//Same as the BVH Scene, but this one construct the BVH with multi-thrading
	//////////////////////////////////////////////////
	class turbo_bvh_scene : public bvh_scene
	{
	public:

		turbo_bvh_scene(Visualizer::Visualizer * visu) :
			bvh_scene(visu)
		{}

	protected:


		virtual void pre_compute_body(int n_max, int ct, int ci, primitive_collection & current_primitives,
			BoundingBox const& current_bb, bvh_type *& current_tree_pos, unsigned int depth)
		{
			assert(current_tree_pos == nullptr);
			assert(current_primitives.size() >= 1);

			//maybe have it depends of the bb 
			double error_margin = 0.001;
			//Math::Vector3f up_layer = (current_bb.max() - current_bb.min()) * error_margin;
			Math::Vector3f up_layer = Math::makeVector(error_margin, error_margin, error_margin);

			BoundingBox tmp_bb = BoundingBox(current_bb.min() - up_layer, current_bb.max() + up_layer);
			if (current_primitives.size() <= n_max)
			{
				current_tree_pos = bvh_type::make_new_leaf(tmp_bb, current_primitives);
			}
			else
			{
				bvh_type * current_node = bvh_type::make_new_node(tmp_bb);
				current_tree_pos = current_node;

				//sparate the current_primitives:

				primitive_collection p_left, p_right;
				BoundingBox bb_left, bb_right;

				find_best_partition_with_sah(ci, ct, current_primitives, current_bb, p_left, p_right, bb_left, bb_right);


				current_primitives.clear();

				if (depth > 3)
				{
					pre_compute_body(n_max, ct, ci, p_left, bb_left, (*current_node)[0], depth + 1);
					pre_compute_body(n_max, ct, ci, p_right, bb_right, (*current_node)[1], depth + 1);
				}
				else
				{
					bvh_type*& left_son = (*current_node)[0];
					bvh_type*& right_son = (*current_node)[1];
					++depth;

					auto left_function = [&](){pre_compute_body(n_max, ct, ci, p_left, bb_left, left_son, depth); };
					auto right_function = [&]() {pre_compute_body(n_max, ct, ci, p_right, bb_right, right_son, depth); };
					
					tbb::parallel_invoke(left_function, right_function);
					
				}
				

			}
		}

	};

}