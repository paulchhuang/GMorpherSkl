/* *************************************************
 * Copyright (2016) : Paul Huang
 * *************************************************/
#ifndef SOLVER_EM_CERES_H_DEFINED
#define SOLVER_EM_CERES_H_DEFINED

// libPatchedMesh
#include <PatchedCloud.h>
// libGMorpher
#include <Solver_ceres.h>
// US
#include <IEnergyTerm_EM_ceres.h>


class Solver_EMCeres
{
	public :
		Solver_EMCeres( bool                                        prediction,
						const SklRgedPatchedCloud::PatchedCloud&           mesh);


		~Solver_EMCeres();


		int solve(const std::list<GMorpher::EPtrCeres>& regularizationTerms,				
	           bool                             probabilistic,
	           int                              maxIter_EM,
	           int                              maxIter_M,
	           int                              maxSubDiv,
	           float                            meanEdge,
	           float                            normThresh,
	           float                            EOutlier,
	           float                            alphaNN,
	           const std::vector<CC3D::float3>& obs_cloud_coords,
	           const std::vector<CC3D::float3>& obs_cloud_normals,
	           std::vector<GMorpher::rigidT>&   RT,
			   std::vector<CC3D::float3>&	    JX,	
	           float&                           sigma);
	

	void setRefMesh( SklRgedPatchedCloud::PatchedCloud&           mesh);

	const GMorpher::SolverCeres& gsolver() const { return *m_GSolver; }

	protected :
		GMorpher::SolverCeres*                 m_GSolver;
	boost::shared_ptr<IEnergyTerm_EMCeres> m_EICP;
};




#endif
