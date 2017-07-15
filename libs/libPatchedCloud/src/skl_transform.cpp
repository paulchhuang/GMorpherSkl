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
// transform skeleton
// ###################################################
// ###################################################

	void transform_ppdicular_ft0( int           numPatches,
								  const float*  RT,	
								  const float3* beta0,
								  float3*       ppdicular_ft )
	{		
		#pragma omp parallel for 
		for(int pi=0;pi<numPatches;++pi) {
			const float* R = RT + 12*pi;
			const float3& t = *(const float3*)(R + 9);			
			ppdicular_ft[pi] = prod3( R, beta0[pi]) + t;		// compare eq 3.19 and 3.13 in the thesis
		}
	}

	void rotate_beta0(	int           numPatches,
						const float*  RT,	
				        const float3* beta0,
				        float3*       beta )
	{		
		#pragma omp parallel for 
		for(int pi=0;pi<numPatches;++pi) {
			const float* R = RT + 12*pi;			
			beta[pi] = prod3( R, beta0[pi]);		// compare eq 3.19 and 3.13 in the thesis
		}
	}

	void transform_JX0_smooth(		const std::vector<int>& jp,
									const std::vector<int>& jp_bound,
									const float*  RT,	
									const float3* dJX0_smooth,
									float3*       JX0_tfed ){

		int numJoints = jp_bound.size() - 1;

		//for (int ji = 0; ji < numJoints; ++ji){
		//	int begin = jp_bound[ji];
		//	int end = jp_bound[ji+1];
		//	int idx = begin;			

		//	const float3* 		dJXi_smooth_itr		= dJX0_smooth + jp_bound[ji];		//begin
		//	const float3* const dJXi_smooth_end		= dJX0_smooth + jp_bound[ji+1];			
		//	
		//	while (dJXi_smooth_itr!=dJXi_smooth_end){
		//		const float* Rj  = RT + 12*jp[idx];
		//		const float3& tj = *(const float3*)(Rj+9);
		//		JX0_tfed[idx++] = prod3(Rj, *dJXi_smooth_itr++) + tj;		
		//	}
		//	assert(idx==end);
		//}

		// Parallel

		#pragma omp parallel for 
		for (int i = 0; i < jp.size(); ++i){
			const float* Rj  = RT + 12*jp[i];		
			const float3& tj = *(const float3*)(Rj+9);
			JX0_tfed[i] = prod3(Rj, dJX0_smooth[i]) + tj;		
		}
		//assert((dJXi_smooth_itr-dJX0_smooth)==jp.size());
	} 

	void rotate_dJX0_smooth(	  const std::vector<int>& jp,
								  const std::vector<int>& jp_bound,
								  const float*  RT,	
								  const float3* dJX0_smooth,
								  float3*       dJX_smooth ){
		
		int numJoints = jp_bound.size() - 1;

		//for (int ji = 0; ji < numJoints; ++ji){
		//	int begin = jp_bound[ji];
		//	int end = jp_bound[ji+1];
		//	int idx = begin;			

		//	const float3* 		dJXi_smooth_itr		= dJX0_smooth + jp_bound[ji];		// begin
		//	const float3* const dJXi_smooth_end		= dJX0_smooth + jp_bound[ji+1];			
		//	
		//	while (dJXi_smooth_itr!=dJXi_smooth_end){
		//		const float* Rj  = RT + 12*jp[idx];				
		//		dJX_smooth[idx++] = prod3(Rj, *dJXi_smooth_itr++);		
		//	}
		//	assert(idx==end);
		//}
		

		// Parallel

		#pragma omp parallel for 
		for (int i = 0; i < jp.size(); ++i){
			const float* Rj  = RT + 12*jp[i];				
			dJX_smooth[i] = prod3(Rj, dJX0_smooth[i]);		
		}
		//assert((dJXi_smooth_itr-dJX0_smooth)==jp.size());
	} 
}