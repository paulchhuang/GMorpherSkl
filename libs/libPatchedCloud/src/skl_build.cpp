/* *************************************************
 * this file describes some functions to build weighting gamma, kappa, trnasform perpendicular foot
 * Copyright (2012) : Paul Huang
 * *************************************************/
#include <PatchedCloud.h>
#include <numeric>
#include <algorithm>
#include <cassert>
#define NUM_JOINTS_SUB 16

namespace SklRgedPatchedCloud 
{

// ###################################################
// ###################################################
// Building skeleton function
// ###################################################
// ###################################################

	void build_joint_coords0(const std::vector<KBNode>&   KBNodes, 
							std::vector<float3>&		 joint_coords0){
		joint_coords0.resize(KBNodes.size());
		for(int ji = 0; ji < KBNodes.size(); ++ji){
			joint_coords0[ji].x = KBNodes[ji].m_x;
			joint_coords0[ji].y = KBNodes[ji].m_y;
			joint_coords0[ji].z = KBNodes[ji].m_z;
		}
	}	

	void build_dJX0_smooth(	const	std::vector<rigidT>& RT0,
							const	std::vector<float3>& JX0, 
							const	std::vector<int>&	jp, 
							const	std::vector<int>&	jp_bound, 
									std::vector<float3>& dJX0_smooth){
	
		int numJoints = jp_bound.size() - 1;
		std::vector<float3> JX0_tmp = JX0;
		
		int offset = 0;
		while ((JX0_tmp.size()%(NUM_JOINTS_SUB-1))!=0){
			std::vector<float3>::iterator JX0_itr = JX0_tmp.begin() + offset;
			JX0_tmp.erase(JX0_itr);				//	after erasing the element, the iterator becomes invalid as well, cannot be reused.
			offset+=(NUM_JOINTS_SUB-1);
		}
		assert(JX0_tmp.size()==numJoints);


		
		dJX0_smooth.resize(jp.size());

		//#pragma omp parallel for 
		for (int ji = 0; ji < numJoints; ++ji){
			int begin = jp_bound[ji];
			int end = jp_bound[ji+1];
			int idx = begin;

			std::vector<float3>::const_iterator dJXi_smooth_begin	= dJX0_smooth.begin() + jp_bound[ji];
			std::vector<float3>::const_iterator dJXi_smooth_end		= dJX0_smooth.begin() + jp_bound[ji+1];
			std::vector<float3>::iterator		dJXi_smooth_itr		= dJX0_smooth.begin() + jp_bound[ji];
			while (dJXi_smooth_itr!=dJXi_smooth_end)
			{
				(*dJXi_smooth_itr) = JX0_tmp[ji] - RT0[jp[idx]].m_t;
				dJXi_smooth_itr++;
				idx++;
			}
			assert(idx==end);
		}
	}

	void build_gamma(	int							numPatches,
						const rigidT*				RT0, 
						const int*					p2j,
						const int*					p2j_child,
						const std::vector<float3>&	j_coords, 
						std::vector<float>&			gamma){
		gamma.resize(numPatches);
		#pragma omp parallel for 
		for(int pi = 0; pi < numPatches; ++pi){
			const float3& ci = RT0[pi].m_t;
			const float3& j = j_coords[p2j[pi]];
			const float3& j_child = j_coords[p2j_child[pi]];

			float norm2Ci_j = (ci.x - j.x)*(ci.x - j.x) + (ci.y - j.y)*(ci.y - j.y) + (ci.z - j.z)*(ci.z - j.z);
			float norm2Ci_j_child = (ci.x - j_child.x)*(ci.x - j_child.x) + (ci.y - j_child.y)*(ci.y - j_child.y) + (ci.z - j_child.z)*(ci.z - j_child.z);
			float norm2j_j_child = (j.x - j_child.x)*(j.x - j_child.x) + (j.y - j_child.y)*(j.y - j_child.y) + (j.z - j_child.z)*(j.z - j_child.z);

			// according to Graz's paper [ECCV12]		
			gamma[pi] = 0.5 - (norm2Ci_j - norm2Ci_j_child)/(2*norm2j_j_child);
		}

	}

	
	void build_kappa(	const std::vector<float>& gamma, 
						const std::vector<int>& pj,
						bool weightflag, std::vector<float>& kappa){
		
		assert(gamma.size()==pj.size());
		kappa.resize(gamma.size());
		
		int numJoints = 2 + (*std::max_element(pj.begin(),pj.end()));

		float sigma = 0.16f;		
		float mu = 0.5f;
		//float kappa_sum = 0;
		std::vector<float> sumWeight(numJoints);
		
		// -------------------------------------
		// 1 - write unnormalised weights and accumulate weights
		for(int pi = 0; pi < kappa.size(); ++pi){
			kappa[pi] = (weightflag) ? exp(- (gamma[pi]-mu)*(gamma[pi]-mu)/(2*sigma*sigma)) : 1;			// using N(0.5, 0.16) as weighting function, s.t. weight at 0 and 1 is 1% of weight at 0.5			
			sumWeight[pj[pi]] += kappa[pi];
		}

		assert(sumWeight[0]==0);

		// -------------------------------------
		// 2 - normalize weights
		#pragma omp parallel for 
		for(int pi = 0; pi < kappa.size(); ++pi)		// then do the normalization such that sum = 1
			kappa[pi] /= sumWeight[pj[pi]];				
	}

	
	
	void build_ppdicular_ft0(	const std::vector<float3>&	j_coords, 
								const int*					p2j, 
								const int*					p2j_child, 
								const std::vector<float>&	gamma, 
								std::vector<float3>&		ppdicular_ft0){

		ppdicular_ft0.resize(gamma.size());
		#pragma omp parallel for 
		for(int pi = 0; pi < ppdicular_ft0.size(); ++pi){			
			const float3& j = j_coords[p2j[pi]];
			const float3& j_child = j_coords[p2j_child[pi]];
			ppdicular_ft0[pi] = gamma[pi] * j + (1 - gamma[pi]) * j_child;
		}

	}

	void build_beta0(	const rigidT*				RT0,
						const std::vector<float3>&	ppdicular_ft0, 
						std::vector<float3>&		beta0){
		
		beta0.resize(ppdicular_ft0.size());
		#pragma omp parallel for 
		for(int pi = 0; pi < ppdicular_ft0.size();++pi){
			const float3& ci = RT0[pi].m_t;
			beta0[pi] = ppdicular_ft0[pi] - ci;
		}
	}
}