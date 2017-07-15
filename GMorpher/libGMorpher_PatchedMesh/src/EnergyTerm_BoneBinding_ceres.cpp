/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#include <EnergyTerm_BoneBinding_ceres.h>
#include <PatchedCloud.h>

#define NUM_JOINTS_SUB 15

namespace SklRgedPatchedCloud{

	
	EnergyTerm_BoneBindingCeres::EnergyTerm_BoneBindingCeres(const PatchedCloud& PC, float EBF) :
	m_pj(PC.patch_to_joint()),
	m_pj_child(PC.patch_to_joint_child()),
	m_numPatches(PC.numPatches()),
	m_EboneWeight(EBF)
	{
		int numPatches = PC.numPatches();
				
		m_beta0				= new float3[numPatches];
		m_w_gamma			= new float[numPatches];
		m_w_kappa			= new float[numPatches];

			
		std::copy(PC.beta0().begin(), PC.beta0().end(), m_beta0 );
		std::copy(PC.w_gamma().begin(), PC.w_gamma().end(), m_w_gamma );
		std::copy(PC.w_kappa().begin(), PC.w_kappa().end(), m_w_kappa );			
		
		#pragma omp parallel for 				
		for(int pi = 0; pi < m_pj.size(); ++pi){		
			m_pj[pi]		-= (1+(m_pj[pi])/(NUM_JOINTS_SUB+1));				// in this function, we don't consider ROOT_T. All elements in p2j and p2j_child minus 1.
			m_pj_child[pi]	-= (1+(m_pj_child[pi])/(NUM_JOINTS_SUB+1));			
		}	

		/*int num_comp = 1 + (*std::max_element(m_pj.begin(), m_pj.end()))/NUM_JOINTS_SUB;
		
		m_componentPatches.resize(num_comp,0);		
		m_p_comp.resize(m_pj.size());
		
		for(std::vector<int>::const_iterator pj_itr = m_pj.begin(); pj_itr != m_pj.end(); ++pj_itr){		
			int cid = (*pj_itr)/NUM_JOINTS_SUB;
			m_componentPatches[cid]++;	
			m_p_comp[pj_itr-m_pj.begin()] = cid;
		}		

		m_componentPatches_accum.resize(num_comp,0);
		for (int ci = 0; ci < num_comp; ++ci)
			m_componentPatches_accum[ci] = std::accumulate(m_componentPatches.begin(),m_componentPatches.begin()+ci+1,0);

		SklRgedPatchedCloud::find_idx_in_jOffsp(m_pj_child, PC.jOffsp(), m_idx_in_jOffsp);*/

	}

	EnergyTerm_BoneBindingCeres::~EnergyTerm_BoneBindingCeres()
	{
		delete[]	m_beta0;
		delete[]	m_w_gamma;
		delete[]	m_w_kappa;
	}

	void EnergyTerm_BoneBindingCeres::addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX)
	{
				
		for(int pi = 0; pi < m_numPatches; ++pi)
		{
			const float* const gamma_ptr = m_w_gamma + pi;
			const float3* const beta0_ptr = m_beta0 + pi;
			
			ceres::CostFunction* cost_function = BetaConstraint::Create(beta0_ptr, gamma_ptr);
			problem.AddResidualBlock(cost_function,
				new ceres::ScaledLoss(NULL, m_EboneWeight*m_w_kappa[pi], ceres::TAKE_OWNERSHIP), 
					RT + 6 * pi, JX + 3 * m_pj[pi], JX + 3 * m_pj_child[pi]);
		}		
	}
	
}
