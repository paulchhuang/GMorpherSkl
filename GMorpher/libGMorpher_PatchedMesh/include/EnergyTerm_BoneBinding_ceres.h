/* *************************************************
* Copyright (2015) : Paul Huang
* *************************************************/
#ifndef ENERGYTERM_BONEBINDING_CERES_H_DEFINED
#define ENERGYTERM_BONEBINDING_CERES_H_DEFINED

#include <IEnergyTerm_ceres.h>
#include "PatchedCloud.h"
namespace SklRgedPatchedCloud
{
	struct BetaConstraint
	{
		BetaConstraint(const float3* const beta0, const float* const gamma) : m_beta0(beta0), m_gamma(double(*gamma))
		{
			m_1_gamma = 1 - m_gamma;
		}

		template <typename T>
		bool operator()(const T* const para_k, const T* const jx, const T* const jx_child, T* residuals) const {

			T beta0[3], pdelta[3];

			beta0[0] = T(m_beta0->x);		beta0[1] = T(m_beta0->y);		beta0[2] = T(m_beta0->z);

			/*current status are computed while computing energy, no need to transform them externally beforehand*/
			ceres::AngleAxisRotatePoint(para_k, beta0, pdelta);

			pdelta[0] += para_k[3];		pdelta[1] += para_k[4];		pdelta[2] += para_k[5];

			residuals[0] = pdelta[0] - (m_gamma*jx[0] + m_1_gamma*jx_child[0]);
			residuals[1] = pdelta[1] - (m_gamma*jx[1] + m_1_gamma*jx_child[1]);
			residuals[2] = pdelta[2] - (m_gamma*jx[2] + m_1_gamma*jx_child[2]);

			return true;
		}

		// Factory to hide the construction of the CostFunction object from the client code.
		static ceres::CostFunction* Create(const float3* const beta0, const float* const gamma)
		{
			return (new ceres::AutoDiffCostFunction<BetaConstraint, 3, 6, 3, 3>(new BetaConstraint(beta0, gamma)));
		}
		const float3* m_beta0;
		double m_gamma;
		double m_1_gamma;
	};

	class EnergyTerm_BoneBindingCeres : public GMorpher::IEnergyTermCeres
	{
	public:
		EnergyTerm_BoneBindingCeres::EnergyTerm_BoneBindingCeres(const PatchedCloud& PC, float EBF);

		virtual ~EnergyTerm_BoneBindingCeres();

		virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX);


	protected:

		float3*					m_beta0;

		// these two set of coefficients remain the same all the time.
		float*					m_w_gamma;				// used for each patch to find a perpendicular foot		
		float*					m_w_kappa;				// used for each patch when summing up bone-binding energy
		float					m_EboneWeight;

		std::vector<int>		m_pj;					// minus 1 here (remove the index of ROOT_T)
		std::vector<int>		m_pj_child;
		std::vector<int>		m_idx_in_jOffsp;
		std::vector<int>		m_componentPatches;
		std::vector<int>		m_componentPatches_accum;
		std::vector<int>		m_p_comp;

		int						m_numPatches;


	};

}

#endif