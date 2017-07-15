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
	void build_w_JX0_smooth(	int numJoint,
								const	std::vector<float3>& dJX0_smooth, 							
								const	std::vector<int>&	jp_bounds, 
										std::vector<float>&	w_JX_smooth){

		std::vector<float>	norm(dJX0_smooth.size());
		std::vector<float>	sumWeight(numJoint);

		#pragma omp parallel for 
		for (int i = 0; i < dJX0_smooth.size(); ++i)
			norm[i] = norm2(dJX0_smooth[i]);				
		
		
		float avg		= std::accumulate(norm.begin(),norm.end(),0.)/((float)dJX0_smooth.size());
		float disVar	= avg*avg;
		
		// -------------------------------------
		// 1 - write unnormalised weights
		#pragma omp parallel for 
		for (int i = 0; i < dJX0_smooth.size(); ++i){
			w_JX_smooth[i] = exp( - (norm[i]*norm[i]) / disVar );			
		}

		// -------------------------------------
		// 2 - accum weights
		//#pragma omp parallel for 
		for (int ji = 0; ji < numJoint; ++ji){			

			std::vector<float>::const_iterator w_begin	= w_JX_smooth.begin() + jp_bounds[ji];
			std::vector<float>::const_iterator w_end	= w_JX_smooth.begin() + jp_bounds[ji+1];
			std::vector<float>::iterator		w_itr	= w_JX_smooth.begin() + jp_bounds[ji];
			while (w_itr!=w_end){
				sumWeight[ji] += (*w_itr);
				w_itr++;				
			}
		}

		// -------------------------------------
		// 3 - normalize weights
		//#pragma omp parallel for 
		for (int ji = 0; ji < numJoint; ++ji){		

			float sumW_ji = sumWeight[ji];
			std::vector<float>::const_iterator w_begin	= w_JX_smooth.begin() + jp_bounds[ji];
			std::vector<float>::const_iterator w_end	= w_JX_smooth.begin() + jp_bounds[ji+1];
			std::vector<float>::iterator		w_itr	= w_JX_smooth.begin() + jp_bounds[ji];
			while (w_itr!=w_end){
				(*w_itr) /= sumW_ji;	
				w_itr++;
			}
		}
	}



	
// ###################################################
// ###################################################
// Inverse skinning thing
// ###################################################
// ###################################################

	void fillBMatrix(	const std::vector<float>&		gamma, 
						const std::vector<float>&		kappa, 
						const std::vector<int>&			pj,
						const std::vector<int>&			pj_child,
							  Eigen::MatrixXf&			B){
		
		B.setZero(B.rows(),B.cols());
		Eigen::Matrix3f I(3,3);		I.setIdentity(3,3);
								  
		std::vector<int> pj_tmp = pj;
		std::vector<int> pj_child_tmp = pj_child;		

		removeROOT_T(pj, pj_tmp);						// in this function, we don't consider ROOT_T. All elements in p2j and p2j_child minus 1.
		removeROOT_T(pj_child, pj_child_tmp);

		#pragma omp parallel for 				
		for(int pi = 0; pi < pj_tmp.size(); ++pi){						
			float w = sqrtf(kappa[pi]);

			B.block(3*pi,3*pj_tmp[pi],3,3)		= w*gamma[pi]*I;
			B.block(3*pi,3*pj_child_tmp[pi],3,3)= w*(1-gamma[pi])*I;
		}

		// debug
		/*std::ofstream f_out("B matrix");	
		f_out<< B;
		f_out.close();*/
	}

	void fillDeltaVec(	const std::vector<float>&	kappa, 
						const std::vector<rigidT>&	RT, 
						const float3*				beta0,
							  Eigen::VectorXf&		Delta){
		
		int numPatches = RT.size();
							
		std::vector<float3> beta0_Rotated(numPatches);

		float* RTf = new float[12*numPatches];
		RigidTFlatten(RT, RTf);

		#pragma omp parallel for 	
		for(int pi=0;pi<numPatches;++pi) {
			const float* R = RTf + 12*pi;			
			beta0_Rotated[pi] = prod3( R, beta0[pi]);		// compare eq 3.19 and 3.13 in the thesis

			float w = sqrtf(kappa[pi]);

			Delta(3*pi+0,0)	= w*(RT[pi].m_t.x + beta0_Rotated[pi].x);
			Delta(3*pi+1,0) = w*(RT[pi].m_t.y + beta0_Rotated[pi].y);
			Delta(3*pi+2,0) = w*(RT[pi].m_t.z + beta0_Rotated[pi].z);
		}

		delete[] RTf;
	}
}