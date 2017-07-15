/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <IEnergyTerm_EM_ceres.h>
#include "kernels_include/kernel_ICPPred.h"
#include <cuda_runtime.h>


EnergyTerm_EMPredCeres::EnergyTerm_EMPredCeres(double                                      alpha,
                                      const SklRgedPatchedCloud::PatchedCloud&           mesh):
IEnergyTerm_EMCeres( alpha, mesh )
{
	// copy topology info
	m_sadj				 = new int[mesh.sadj().size()];
	m_sadj_bounds		 = new int[m_numPatches+1];
	m_smooth_bounds_host = new int[m_numPatches+1];
	
	std::copy( mesh.sadj().begin(), mesh.sadj().end(), m_sadj);
	std::copy( mesh.sadj_bounds().begin(), mesh.sadj_bounds().end(), m_sadj_bounds);
	std::copy( mesh.smooth_bounds().begin(), mesh.smooth_bounds().end(), m_smooth_bounds_host);

	// copy the dXN0
	m_numVertices_smooth = mesh.numVertices_smooth();
	m_dXN0_smooth    = new float3[2*m_numVertices_smooth];
	m_XN_smooth_host = new float3[2*m_numVertices_smooth];

	reinit(mesh); //this will have been called 

	cudaMalloc( &m_XN_smooth_device, sizeof(float3)*2*m_numVertices_smooth );
	cudaMalloc( &m_smooth_bounds_device, sizeof(int)*(m_numPatches+1) );
	cudaMemcpy( m_smooth_bounds_device, &(mesh.smooth_bounds()[0]), sizeof(int)*(m_numPatches+1), cudaMemcpyHostToDevice );
}



void EnergyTerm_EMPredCeres::reinit(const SklRgedPatchedCloud::PatchedCloud& mesh)
{
	// first reinitialise the base mesh
	IEnergyTerm_EMCeres::reinit(mesh);

	// second reinitialize the rest
	const float3*  dX_itr = &(mesh.dX0_smooth()[0]);
	float3*       dXN_itr = m_dXN0_smooth;
	for(int pi=0;pi<m_numPatches;++pi) {
		int numVertices_i   = m_pv_bounds_host[pi+1] - m_pv_bounds_host[pi];
		int numNeighbours_i = m_sadj_bounds[pi+1] - m_sadj_bounds[pi];
		const std::vector<float3>::const_iterator N_begin = mesh.N0().begin() + m_pv_bounds_host[pi];
		for(int ni=0;ni<numNeighbours_i+1;++ni){
			std::vector<float3>::const_iterator N_itr = N_begin;
			for(int vi=0;vi<numVertices_i;++vi) {
				dXN_itr[0] = *dX_itr++;
				dXN_itr[1] = *N_itr++;
				dXN_itr+=2;
			}
		}
	}
}

EnergyTerm_EMPredCeres::~EnergyTerm_EMPredCeres()
{
	delete[] m_sadj;
	delete[] m_sadj_bounds;
	delete[] m_smooth_bounds_host;
	delete[] m_dXN0_smooth;
	delete[] m_XN_smooth_host;

	cudaFree( m_XN_smooth_device );
	cudaFree( m_smooth_bounds_device );	
}






void EnergyTerm_EMPredCeres::EStep(const float                     sigma,
                               const float                     normThresh,
                               const float                     EOutlier,
                               const float*                    RT )
{
	// transform cloud and upload to GPU
	SklRgedPatchedCloud:: transform_Cloudsmooth( m_numPatches, RT, m_smooth_bounds_host, m_sadj, m_sadj_bounds, m_pv_bounds_host, m_dXN0_smooth, m_XN_smooth_host );
	cudaMemcpy( m_XN_smooth_device, m_XN_smooth_host, sizeof(float3)*2*m_numVertices_smooth, cudaMemcpyHostToDevice );

//	cudaError_t status = cudaGetLastError();
//	if(status != cudaSuccess)std::cout<<__FUNCTION__<<":"<<__LINE__<<" "<<cudaGetErrorString(status)<<"\n";


	// run the heavy find nearest neighbor kernel

	runKernel_ICPPred_EStep( m_numPatches,
	                         m_numObs,
	                         sigma,
	                         normThresh,
	                         EOutlier,
	                         m_XNObs_device,
	                         m_XN_smooth_device,
	                         m_smooth_bounds_device,
	                         m_weights_device,
	                         m_indices_device);


	// get the stuff back
	cudaMemcpy( (void*)m_weights_host, m_weights_device, sizeof(float)*m_numObs*(m_numPatches+1), cudaMemcpyDeviceToHost );
	cudaMemcpy( (void*)m_indices_host, m_indices_device, sizeof(int)*m_numObs*(m_numPatches+1), cudaMemcpyDeviceToHost );


	// run the return gymnastic on indices 
	// we get min_ptr - beg_ptr; so we need to 
	// modulo numVertices_i * 2
	// add jv_bounds * 2
	#pragma omp parallel for 
	for(int pi=0;pi<m_numPatches;++pi){
		int numVertices_ix2 = (m_pv_bounds_host[pi+1] - m_pv_bounds_host[pi])*2;
		if( numVertices_ix2 == 0 ) numVertices_ix2 = 1; // in case the freaking bone has 0 vertices associated
		int offset          = m_pv_bounds_host[pi]*2;
		int* i_ptr = m_indices_host + pi;
		for(int vi=0;vi<m_numObs;++vi) {
			*i_ptr = (*i_ptr % numVertices_ix2) + offset;
			i_ptr += m_numPatches + 1; 
		}
	}
}
