/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef ENERGYTERM_RIGIDITY_CEDRIC_CERES_H_DEFINED
#define ENERGYTERM_RIGIDITY_CEDRIC_CERES_H_DEFINED

#include <IEnergyTerm_ceres.h>
#include "PatchedCloud.h"

namespace SklRgedPatchedCloud
{

	struct RigidityConstraintCEDRIC
	{
		RigidityConstraintCEDRIC(const float3* const dX_k, const float3* const dX_l) : m_dX_k(dX_k), m_dX_l(dX_l)	{}

		template <typename T>
		bool operator()(const T* const para_k, const T* const para_l, T* residuals) const {

			T pk0[3], pl0[3];
			T pk[3], pl[3];

			pk0[0] = T(m_dX_k->x);		pk0[1] = T(m_dX_k->y);		pk0[2] = T(m_dX_k->z);
			pl0[0] = T(m_dX_l->x);		pl0[1] = T(m_dX_l->y);		pl0[2] = T(m_dX_l->z);

			/*current status are computed while computing energy, no need to transform them externally beforehand*/
			ceres::AngleAxisRotatePoint(para_k, pk0, pk);
			ceres::AngleAxisRotatePoint(para_l, pl0, pl);

			pk[0] += para_k[3];		pl[0] += para_l[3];
			pk[1] += para_k[4];		pl[1] += para_l[4];
			pk[2] += para_k[5];		pl[2] += para_l[5];

			residuals[0] = pk[0] - pl[0];
			residuals[1] = pk[1] - pl[1];
			residuals[2] = pk[2] - pl[2];

			return true;
		}

		// Factory to hide the construction of the CostFunction object from the client code.
		static ceres::CostFunction* Create(const float3* const pPoint_k, const float3* const pPoint_l)
		{
			return (new ceres::AutoDiffCostFunction<RigidityConstraintCEDRIC, 3, 6, 6>(new RigidityConstraintCEDRIC(pPoint_k, pPoint_l)));
		}
		const float3* m_dX_k;
		const float3* m_dX_l;
	};


	class EnergyTerm_Rigidity_CEDRIC_ceres : public GMorpher::IEnergyTermCeres
	{
		public :
			EnergyTerm_Rigidity_CEDRIC_ceres(const PatchedCloud&  PC);

			virtual ~EnergyTerm_Rigidity_CEDRIC_ceres();

	
			virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX);

		protected :
		float3* m_dX0_smooth;						
		float*  m_w_smooth;

		std::vector<int> m_sadj;
		std::vector<int> m_sadj_bounds;
		std::vector<int> m_pv_bounds;
		std::vector<int> m_smooth_bounds;		
		std::vector<int> m_patch_comp;	
	};

} // end namespace PatchedCloud


#endif
