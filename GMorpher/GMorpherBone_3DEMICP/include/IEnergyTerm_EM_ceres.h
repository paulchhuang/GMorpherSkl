/* *************************************************
 * Copyright (2016) : Paul Huang
 * *************************************************/
#ifndef IENERGYTERM_EM_CERES_H_DEFINED
#define IENERGYTERM_EM_CERES_H_DEFINED

// libGMorpher
#include <IEnergyTerm_ceres.h>
// libPatchedMesh
#include <PatchedCloud.h>

class IEnergyTerm_EMCeres : public GMorpher::IEnergyTermCeres
{
	public :
	typedef  SklRgedPatchedCloud::float3 float3;
	IEnergyTerm_EMCeres(double                                      alpha,
	                const SklRgedPatchedCloud::PatchedCloud&           mesh);


	virtual ~IEnergyTerm_EMCeres();


	/// this simply resizes the target cloud and the weights/targets matrices
	void setFrameData( const std::vector<CC3D::float3>& obs_cloud_coords,
	                   const std::vector<CC3D::float3>& obs_cloud_normals );

	inline void setAlpha( double alpha ) { m_alpha = alpha; }

	virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX);

	/*virtual double computeEnergy( const float* RT, const float* JX);
	virtual void   addToGTG_GTb( const float* RT, const float* JX, GMorpher::BigMatrix& BM);*/
	double reEvaluateSigma( const float* RT )const;


	/// this reevaluates the weights and index matrices
	virtual void EStep( const float                     sigma,
	                    const float                     normThresh,
	                    const float                     EOutlier,
	                    const float*                    RT ) = 0;

	/// this modifies the weights and index matrices so only one patch per 
	/// vertex remains, simulating the behavior of a normal ICP algorithm.
	virtual void EStep_makeDeterministicAssignment();

	virtual void reinit( const SklRgedPatchedCloud::PatchedCloud& mesh);

	void setP2J(const std::vector<int>& p2j){ m_p2j = p2j; }

	protected :
	double m_alpha;
	// the capacity of the smooth vector
	int m_numVertices;
	// the number of patches expected
	int m_numPatches;


	// both have to be under maxNumTargets
	int               m_numObs;	
	float3*           m_XNObs_host;
	void*             m_XNObs_device;
	void*             m_weights_device;  // floats (numPatches+1)*numTargets	
	float*            m_weights_host;    // floats (numPatches+1)*numTargets
	void*             m_indices_device;  // int (numPatches+1)*numTargets
	int*              m_indices_host;    // int (numPatches+1)*numTargets

	float3*  m_XN0;
	float3*  m_dXN0_host;				// 2*numVertices
	float3*  m_XN_host;					// 2*numVertices
	int*     m_pv_bounds_host;         //int (numPatches+1)

	std::vector<int> m_patch_component;	
	std::vector<int> m_p2j;
};


//class EnergyTerm_EMCeres : public IEnergyTerm_EMCeres
//{
//	public :
//
//		EnergyTerm_EMCeres(double                                      alpha,
//	              const SklRgedPatchedCloud::PatchedCloud&           mesh);
//		virtual ~EnergyTerm_EMCeres();
//
//
//	virtual void EStep( const float                     sigma,
//	                    const float                     normThresh,
//	                    const float                     EOutlier,
//	                    const float*                    RT );
//
//	protected :
//	void*   m_XN_device;        // float 6*numVertices
//	void*   m_pv_bounds_device; //int (numPatches+1)
//};




class EnergyTerm_EMPredCeres : public IEnergyTerm_EMCeres
{
	public :

		EnergyTerm_EMPredCeres(double                                      alpha,
	                  const SklRgedPatchedCloud::PatchedCloud&           mesh);

		EnergyTerm_EMPredCeres(double                                      alpha,
                       const SklRgedPatchedCloud::PatchedCloud&    mesh, 
						const std::vector<int>&					  p2j);

		virtual ~EnergyTerm_EMPredCeres();


	virtual void EStep( const float                     sigma,
	                    const float                     normThresh,
	                    const float                     EOutlier,
	                    const float*                    RT );

	virtual void reinit( const SklRgedPatchedCloud::PatchedCloud& mesh);

	protected :
	int m_numVertices_smooth;

	int*    m_sadj;
	int*    m_sadj_bounds;				//int (numPatches+1)

	int*	m_smooth_bounds_host;		//int (numPatches+1)

	float3* m_dXN0_smooth;
	float3* m_XN_smooth_host;
	
	void*   m_XN_smooth_device;         // float 6*numVertices_smooth
	void*   m_smooth_bounds_device;		//int (numPatches+1)
};






#endif
