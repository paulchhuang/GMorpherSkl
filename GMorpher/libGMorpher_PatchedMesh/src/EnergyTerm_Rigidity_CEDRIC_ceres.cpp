/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <EnergyTerm_Rigidity_CEDRIC_ceres.h>
#include <GraphUtils.h>

namespace SklRgedPatchedCloud
{
	
EnergyTerm_Rigidity_CEDRIC_ceres::
EnergyTerm_Rigidity_CEDRIC_ceres(const PatchedCloud&        PC):
m_sadj         (PC.sadj()),
m_sadj_bounds  (PC.sadj_bounds()),
m_pv_bounds    (PC.pv_bounds()),
m_smooth_bounds(PC.smooth_bounds())
{

	int numPatches = m_pv_bounds.size() - 1;
	int numVertices = m_pv_bounds[numPatches];
	int numVertices_smooth = m_smooth_bounds[numPatches];

	m_dX0_smooth = new float3[numVertices_smooth];		
	m_w_smooth   = new float[numVertices_smooth];
	std::copy( PC.dX0_smooth().begin(), PC.dX0_smooth().end(), m_dX0_smooth );
	std::copy( PC.w_smooth().begin(), PC.w_smooth().end(), m_w_smooth );
	
	
	//int numComponents = GMorpher::getConnectedComponents( m_sadj, m_sadj_bounds, m_patch_comp);
	
}

EnergyTerm_Rigidity_CEDRIC_ceres::~EnergyTerm_Rigidity_CEDRIC_ceres()
{
	delete[] m_dX0_smooth;	
	delete[] m_w_smooth;
}

void EnergyTerm_Rigidity_CEDRIC_ceres::addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX)
{
	int numPatches = m_pv_bounds.size() - 1;	
	
//#pragma omp parallel for
	for (int pi = 0; pi < numPatches; ++pi)
	{
		int numVertices_i = m_pv_bounds[pi + 1] - m_pv_bounds[pi];
		
		const std::vector<int>::const_iterator nBegin = m_sadj.begin()	+ m_sadj_bounds[pi];
		const std::vector<int>::const_iterator nEnd = m_sadj.begin()	+ m_sadj_bounds[pi + 1];

		const float3* const dXii_begin = m_dX0_smooth + m_smooth_bounds[pi];
		const float*  const wii_begin = m_w_smooth + m_smooth_bounds[pi];
		const float3*		dXii_end = dXii_begin + numVertices_i;
		const float3*       dXji_itr = dXii_end;
		const float*        wji_itr = wii_begin + numVertices_i;

		for (std::vector<int>::const_iterator pj_itr = nBegin; pj_itr != nEnd;	++pj_itr) 
		{
			int pj = *pj_itr;
			
			const float*        wii_itr = wii_begin;
			const float3*       dXii_itr = dXii_begin;						

			while (dXii_itr != dXii_end)
			{
				
				ceres::CostFunction* cost_function = RigidityConstraintCEDRIC::Create(dXii_itr, dXji_itr);
				problem.AddResidualBlock(cost_function,
					new ceres::ScaledLoss(NULL, (*wji_itr + *wii_itr), ceres::TAKE_OWNERSHIP),/* squared loss */
					RT + 6 * pi, RT + 6 * pj);
				
				dXii_itr++;
				dXji_itr++;
				wji_itr++;
				wii_itr++;
			}
		}
	}
	
}

} // end namespace PatchedCloud
