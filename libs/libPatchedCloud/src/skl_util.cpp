/* *************************************************
 * this file describes some functions to build weighting gamma, kappa, trnasform perpendicular foot
 * Copyright (2012) : Paul Huang
 * *************************************************/
#include <PatchedCloud.h>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <cassert>
#define NUM_JOINTS_SUB 16

namespace SklRgedPatchedCloud 
{	

// ###################################################
// ###################################################
// AUX function
// ###################################################
// ###################################################

	int getChildId(const int ji, const std::vector<KBJoint>& KBTree){
		// this function works for all the joints except ROOT_ROTATE, which has 5 children NECK, RHIP, LHIP, RSHOULDER, LSHOULDER
		for(int joint = 0; joint < KBTree.size(); ++joint){
			if(KBTree[joint].m_parentId==ji){	
				return joint;}
		}
		std::cout<<"there are some problems when find child for joint "<< ji <<std::endl;		
		return EXIT_FAILURE;
	}
	
	int getChildIdKNN(const float3& ci,  const std::vector<KBJoint>& KBTree, const std::vector<float3>&	j_coords){
		// Entering this function indicates that we're finding the closet child of ROOT_ROTATE for the current patch with mean ci
		int closest_child_id;
		float min = FLT_MAX ;
		float3 j_coords_tmp;
		for(int ji = 0; ji < KBTree.size(); ++ji){
			if((KBTree[ji].m_parentId%NUM_JOINTS_SUB)==1){
				j_coords_tmp = j_coords[ji];
				float norm2_tmp = (ci.x - j_coords_tmp.x)*(ci.x - j_coords_tmp.x) + (ci.y - j_coords_tmp.y)*(ci.y - j_coords_tmp.y) + (ci.z - j_coords_tmp.z)*(ci.z - j_coords_tmp.z);
				if (norm2_tmp < min){
					min = norm2_tmp;
					closest_child_id = ji;
				}
			}
		}
		return closest_child_id;
	}

	inline void removeROOT_T(const std::vector<int>& before, std::vector<int>& after){
		
		int size = before.size();		
		assert(after.size()==size);

		#pragma omp parallel for 				
		for(int pi = 0; pi < size; ++pi)		
			after[pi] -= (1+(after[pi])/(NUM_JOINTS_SUB));				// in this function, we don't consider ROOT_T. All elements in p2j and p2j_child minus 1.					
	}


	void compute_patch_to_joint( const SklRgedPatchedCloud::Patching&	patches,
								 const RiggedMesh&				mesh,
								 std::vector<int>&						patch_to_joint )
	{
		// resize
		patch_to_joint.resize(patches.numPatches() );

		// do the computation
		for(int pi=0; pi < patches.numPatches(); ++pi ) {
			// allocate a small buffer to count joint occurences on the patches
			std::vector<int> accum( mesh.numJoints(), 0 );
			std::vector<int>::const_iterator vi_itr = patches.pv().begin()
															  + patches.pv_boundaries()[pi];
			std::vector<int>::const_iterator vi_end = patches.pv().begin()
															+ patches.pv_boundaries()[pi+1];
			// count jointId occurences on all of the vertices in the patch
			while(vi_itr != vi_end ) {
				// these vi_itr indices are indices from the riggedMesh in its 
				// coords form (organized by joint)
				std::vector<int>::const_iterator f_itr =  std::upper_bound(
						 mesh.jv_bounds().begin(), mesh.jv_bounds().end(), *vi_itr);
				// upper_bound returns the index of the first that compares greater than value
				// so that if *vi_itr = 3 and mesh.jv_bounds() = [0,2,4]
				// f_itr - mesh.jv_bounds().begin() will be 2
				// so we still need to subtract 1
				int jointId = f_itr - mesh.jv_bounds().begin() - 1;
				// accumulate for this joint
				accum[jointId]++;
				vi_itr++;
			}

			// select the most frequent
			std::vector<int>::const_iterator f_itr =
									   std::max_element(accum.begin(), accum.end());
			int majJoint = f_itr - accum.begin();
			if( accum[majJoint] == 0 ) throw std::runtime_error("empty patch"); // should only happen on an empty patch
			patch_to_joint[pi] = majJoint;
		}
	}

	void compute_patch_to_joint_child(	int							numPatches,
										const rigidT*				RT0,
										const std::vector<KBJoint>& KBTree, 
										const std::vector<float3>&	j_coords, 
										const int*					p2j, 									
										std::vector<int>&			p2j_child){
		p2j_child.resize(numPatches);	

		for(int pi = 0; pi < numPatches; ++pi){
			int child_id;
			float3 ci = RT0[pi].m_t;
			if((p2j[pi]%NUM_JOINTS_SUB)==1){
				child_id = getChildIdKNN(ci,KBTree,j_coords);}
			else{
				child_id = getChildId(p2j[pi],KBTree);}
			p2j_child[pi] = child_id;	
		}
	}

	void compute_jointOffspring(const std::vector<KBJoint>& KBTree, std::vector<int>& jOffsp, std::vector<int>& jOffsp_bounds){
		// jOffsp_bounds, jOffsp start with ROOT_R (1), not ROOT_T (0)
		
		int numComps = KBTree.size()/NUM_JOINTS_SUB;
		
		jOffsp.reserve(NUM_JOINTS_SUB - 1);
		jOffsp_bounds.resize(KBTree.size()-numComps+1);

		jOffsp_bounds[0] = 0;
		for(int curr_ji = 1; curr_ji < NUM_JOINTS_SUB; ++curr_ji ){		// notice that here we start with ji = 1 (ROOT_R), not ji = 0
			for(int ji = 1; ji < KBTree.size(); ++ji){
				if(KBTree[ji].m_parentId==curr_ji){
					jOffsp.push_back(ji);
					jOffsp_bounds[curr_ji]++;
				}
			}
			jOffsp_bounds[curr_ji] += jOffsp_bounds[curr_ji-1];
		}		

		jOffsp.resize(numComps*(NUM_JOINTS_SUB - 2));
		for (int ci = 1; ci < numComps; ++ci){
			int start	= ci*(NUM_JOINTS_SUB - 2);
			int end		= (ci+1)*(NUM_JOINTS_SUB - 2);	
			for (int i = start; i < end; ++i)			
				jOffsp[i] = jOffsp[i-(NUM_JOINTS_SUB - 2)] + NUM_JOINTS_SUB;			
		}

		for (int ci = 1; ci < numComps; ++ci){
			int start	= ci*(NUM_JOINTS_SUB - 1)+1;
			int end		= (ci+1)*(NUM_JOINTS_SUB - 1)+1;	
			for (int i = start; i < end; ++i)			
				jOffsp_bounds[i] = jOffsp_bounds[i-(NUM_JOINTS_SUB - 1)] + (NUM_JOINTS_SUB - 2);
		}

	}

	void find_idx_in_jOffsp(const std::vector<int>& p2j_child, const std::vector<int>& jOffsp, std::vector<int>& idx_in_jOff){

		/*find the idx of p2j_child in jOffsp, it is useful when accessing r(1-r) matrix in add2GTG_Gb function */		
		int numPatches = p2j_child.size();
		std::vector<int> jOffsp_shifted = jOffsp;
		removeROOT_T(jOffsp, jOffsp_shifted);

		//for(int i=0; i < jOffsp_shifted.size(); ++i)	jOffsp_shifted[i] -= (1+(jOffsp_shifted[i])/NUM_JOINTS_SUB);

		assert(*(std::min_element(p2j_child.begin(), p2j_child.end()))==1);
		assert(*(std::min_element(jOffsp_shifted.begin(), jOffsp_shifted.end()))==1);
		idx_in_jOff.resize(numPatches);

		#pragma omp parallel for 
		for(int pi=0; pi < numPatches; ++pi){
			const std::vector<int>::const_iterator f_itr = std::find(jOffsp_shifted.begin(), jOffsp_shifted.end(), p2j_child[pi]);
			idx_in_jOff[pi] = f_itr - jOffsp_shifted.begin();
		}

	}

	void build_joint_to_patch(	int numJoint, 
								const std::vector<int>&	pj, 
								const std::vector<int>&	pj_child, 
									  std::vector<int>&	jp, 
									  std::vector<int>& jp_bound){
		jp.reserve(pj.size());
		jp_bound.resize(numJoint+1);

		std::vector<int> pj_tmp = pj;		
		std::vector<int> pj_child_tmp = pj_child;

		//#pragma omp parallel for 				
		//for(int pi = 0; pi < pj_tmp.size(); ++pi){		
		//	pj_tmp[pi]		-= (1+(pj_tmp[pi])/(NUM_JOINTS_SUB));				// in this function, we don't consider ROOT_T. All elements in p2j and p2j_child minus 1.
		//	pj_child_tmp[pi]-= (1+(pj_child_tmp[pi])/(NUM_JOINTS_SUB));
		//}

		removeROOT_T(pj, pj_tmp);
		removeROOT_T(pj_child, pj_child_tmp);
		
		jp_bound[0] = 0;
		for (int ji = 0; ji < numJoint; ++ji){
			jp_bound[ji+1] = jp_bound[ji];

			// -------------------------------------
			// 1 - patches attached to ji. 
			std::vector<int>::iterator pj_itr = pj_tmp.begin();
			do{
				pj_itr = std::find(pj_itr, pj_tmp.end(), ji);
				if(pj_itr!=pj_tmp.end()){
					jp.push_back(pj_itr-pj_tmp.begin());
					pj_itr ++;
					jp_bound[ji+1]++;
				}
			}while (pj_itr!=pj_tmp.end());	


			// -------------------------------------
			// 2 - patches attached to ji's parent, i.e. p2j_child = ji
			pj_itr = pj_child_tmp.begin();
			do{
				pj_itr = std::find(pj_itr, pj_child_tmp.end(), ji);
				if(pj_itr!=pj_child_tmp.end()){
					jp.push_back(pj_itr-pj_child_tmp.begin());
					pj_itr ++;
					jp_bound[ji+1]++;
				}
			}while (pj_itr!=pj_child_tmp.end());	
		}

		//assert(jp.size()==jp.capacity());
		jp.resize(jp.size());
	}

// ###################################################
// ###################################################
// Output function (for check)
// ###################################################
// ###################################################
	void outputPatchCenters(int numPatches, const char* filename, const rigidT* RT){
		
		std::ofstream f_out(filename);				
		for (int pi=0;pi<numPatches;++pi){
			const float3&   ci = RT[pi].m_t;			
			f_out<<ci.x<<" "<<ci.y<<" "<<ci.z<<"\n";
		}
		

		f_out.close();
	}

	void outputJoints(int numJoints, const char* filename, const std::vector<float3>&	j_coords){
		std::ofstream f_out(filename);				
		if(!f_out.is_open())throw(std::ios_base::failure(filename));
		
		for (int ji=0;ji<numJoints;++ji){
			const float3&  ci = j_coords[ji];			
			f_out<<ci.x<<" "<<ci.y<<" "<<ci.z<<"\n";
		}		

		f_out.close();
	}

	void outputGamma(int numPatches, const char* filename, const std::vector<float>& gamma){
		
		std::ofstream f_out(filename);				
		for (int pi=0;pi<numPatches;++pi){		
			f_out<< gamma[pi] <<"\n";
		}
		
		f_out.close();
	}

	void outputPPfts(int numPatches, const char* filename, const std::vector<float3>&  ppdicular_ft0){
		
		std::ofstream f_out(filename);				
		for (int pi=0;pi<numPatches;++pi){		
			const float3&  ci = ppdicular_ft0[pi];	
			f_out<<ci.x<<" "<<ci.y<<" "<<ci.z<<"\n";
		}
		
		f_out.close();
	}

	void outputBeta(int numPatches, const char* filename, const std::vector<float3>&  beta0){
		
		std::ofstream f_out(filename);				
		for (int pi=0;pi<numPatches;++pi){		
			const float3&  ci = beta0[pi];	
			f_out<<ci.x<<" "<<ci.y<<" "<<ci.z<<"\n";
		}
		
		f_out.close();
	}
}