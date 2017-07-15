/* *************************************************
 * Copyright (2016) : Paul Huang
 * *************************************************/
#include <IEnergyTerm_EM_ceres.h>
#include <EnergyTerm_3DConstraint_ceres.h>
#include "ICP_function.h"
#include <GraphUtils.h>
#include <cuda_runtime.h>


IEnergyTerm_EMCeres::IEnergyTerm_EMCeres(double                                      alpha,
                                const SklRgedPatchedCloud::PatchedCloud&           mesh)
{
	m_alpha       = alpha;
	m_numPatches  = mesh.numPatches();
	m_numVertices = mesh.numVertices();

	m_XN0        = new float3[2*m_numVertices]; 
	m_dXN0_host  = new float3[2*m_numVertices]; 
	m_XN_host    = new float3[2*m_numVertices];   
	m_pv_bounds_host = new int[m_numPatches+1];

	// copy bounds
	std::copy(mesh.pv_bounds().begin(), mesh.pv_bounds().end(), m_pv_bounds_host);
	// copy coords
	//reinit(mesh); this must be done by the derived class

	// all variables depending on obs
	m_numObs            = 0;
	m_XNObs_host        = 0;
	m_XNObs_device      = 0;
	m_weights_device    = 0;
	m_weights_host      = 0; 
	m_indices_device    = 0; 
	m_indices_host      = 0; 
	
	int numComponents	= GMorpher::getConnectedComponents( mesh.sadj(), mesh.sadj_bounds(), m_patch_component);
}

void IEnergyTerm_EMCeres::addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX)
{
	const float3* xnobs_ptr = m_XNObs_host;
	const float*  w_ptr = m_weights_host;
	const int*    i_ptr = m_indices_host;
	
	for (int ti = 0; ti < m_numObs; ++ti){

		xnobs_ptr = m_XNObs_host + ti * 2;
		i_ptr = m_indices_host + ti*(m_numPatches + 1);
		w_ptr = m_weights_host + ti*(m_numPatches + 1);

		for (int pi = 0; pi < m_numPatches; ++pi){

			int vi = i_ptr[pi];
			if (w_ptr[pi] >= 1e-7)
			{
				ceres::CostFunction* cost_function = GMorpher::Point3DConstraint::Create(m_dXN0_host + vi, xnobs_ptr);
				problem.AddResidualBlock(cost_function,
					new ceres::ScaledLoss(NULL, m_alpha*w_ptr[pi], ceres::TAKE_OWNERSHIP),/* squared loss */
					RT + 6 * pi);
			}
			
		}
	}
}

void IEnergyTerm_EMCeres::reinit(const SklRgedPatchedCloud::PatchedCloud& mesh)
{
	assert( mesh.numVertices() == m_numVertices);
	// copy XN0 and dXN0
	for(int vi=0;vi<m_numVertices;++vi) {
		m_XN0[2*vi]   = mesh.X0()[vi];
		m_XN0[2*vi+1] = mesh.N0()[vi];
	}
	SklRgedPatchedCloud::build_dXN0( m_numPatches, &(mesh.RT0()[0]), m_pv_bounds_host, 
	                                                       m_XN0, m_dXN0_host );
}



IEnergyTerm_EMCeres::~IEnergyTerm_EMCeres()
{
	delete[] m_XN0;
	delete[] m_dXN0_host;
	delete[] m_XN_host;
	delete[] m_pv_bounds_host;

	if( m_numObs != 0 ) {
		delete[] m_XNObs_host;
		delete[] m_weights_host;
		delete[] m_indices_host;
		cudaFree(m_XNObs_device);
		cudaFree(m_weights_device);
		cudaFree(m_indices_device);
	}
}


void IEnergyTerm_EMCeres::setFrameData(const std::vector<CC3D::float3>& obs_cloud_coords,
                                   const std::vector<CC3D::float3>& obs_cloud_normals )
{
	if( m_numObs != 0 ) {
		delete[] m_XNObs_host;
		delete[] m_weights_host;
		delete[] m_indices_host;
		cudaFree(m_XNObs_device);
		cudaFree(m_weights_device);
		cudaFree(m_indices_device);
	}

	m_numObs         = obs_cloud_coords.size();
	// host alloc
	m_XNObs_host     = new float3[2*m_numObs];
	m_weights_host   = new  float[m_numObs*(m_numPatches+1)];
	m_indices_host   = new    int[m_numObs*(m_numPatches+1)];
	// device alloc
	cudaMalloc(&m_XNObs_device,   sizeof(float3)*2*m_numObs );
	cudaMalloc(&m_weights_device, sizeof(float)*m_numObs*(m_numPatches+1) );
	cudaMalloc(&m_indices_device, sizeof(int)*m_numObs*(m_numPatches+1) );	

	// copy the data to the host
	CC3D::float3* v = (CC3D::float3*)m_XNObs_host;
	
	#pragma omp parallel for 
	for(int vi=0; vi<m_numObs; ++vi) { // careful here we dont use the padded numTargets
		/*v[0] = obs_cloud_coords[vi];
		v[1] = obs_cloud_normals[vi];
		v+=2;*/

		v[2*vi+0] = obs_cloud_coords[vi];
		v[2*vi+1] = obs_cloud_normals[vi];
	}

	// copy the data to the device
	cudaMemcpy( m_XNObs_device, (const void*) m_XNObs_host, 
	                        sizeof(float3)*2*m_numObs, cudaMemcpyHostToDevice );
}


void IEnergyTerm_EMCeres::EStep_makeDeterministicAssignment()
{
	float*  w_ptr   = m_weights_host;

	for (int ti = 0; ti < m_numObs; ++ti) 
	{
		float* w_end = w_ptr + (m_numPatches); // HACK... could be +1 if we authorized outliers
		float* f_itr = std::max_element(w_ptr, w_end);
		float val = *f_itr;
		std::fill(w_ptr, w_end, 0.); // set the whole line to 0
		if( val != 0 ) *f_itr = 1.0; // set the best guy to 1 only if we had found a neighbor (if not the kernel puts 0 in w)
		w_ptr += (m_numPatches+1);
	}

//	exit(0);
}

double IEnergyTerm_EMCeres::reEvaluateSigma(const float* RT) const
{	
	
	
	// create the normal cloud from the current state
	SklRgedPatchedCloud::transform_Cloud( m_numPatches, RT, m_pv_bounds_host, 
	                                                    m_dXN0_host, m_XN_host);

	double E = 0.;
	double sumW = 0.;
	const float3* xno_ptr = m_XNObs_host;
	const float*  w_ptr   = m_weights_host;
	const int*    i_ptr   = m_indices_host;

	/*for (int ti = 0; ti < m_numObs; ++ti) 
	{
		const float3& xo = xno_ptr[0];
		const float3& no = xno_ptr[1];
		for(int pi = 0; pi < m_numPatches; ++pi) 
		{
			int vi = i_ptr[pi];
			const float3& x = m_XN_host[vi];
			const float3& n = m_XN_host[vi+1];
			const float   w = w_ptr[pi];
			accumEnergy( w, xo, no, x, n, E );
			sumW += w;
		}
		xno_ptr += 2;
		i_ptr += (m_numPatches+1);
		w_ptr += (m_numPatches+1);
	}
	
	std::cout<< "E: " << E << std::endl;
	std::cout<< "sumW: " << sumW << std::endl;*/


	std::vector<double> EObs(m_numObs);
	std::vector<double> WObs(m_numObs);

	#pragma omp parallel for
	for (int ti = 0; ti < m_numObs; ++ti) 
	{
		xno_ptr = m_XNObs_host + ti*2;
		i_ptr	= m_indices_host + ti*(m_numPatches+1);
		w_ptr	= m_weights_host + ti*(m_numPatches+1);

		const float3& xo = xno_ptr[0];
		const float3& no = xno_ptr[1];
		for(int pi = 0; pi < m_numPatches; ++pi) 
		{
			int vi = i_ptr[pi];
			const float3& x = m_XN_host[vi];
			const float3& n = m_XN_host[vi+1];
			float   w = w_ptr[pi];

			accumEnergy( w, xo, no, x, n, EObs[ti] );
			WObs[ti] += w;
		}
	}	


	//return 0.5* sqrt( std::accumulate(EObs.begin(), EObs.end(),0.0) / (3.0*std::accumulate(WObs.begin(), WObs.end(),0.0)) );
	return sqrt( std::accumulate(EObs.begin(), EObs.end(),0.0) / (3.0*std::accumulate(WObs.begin(), WObs.end(),0.0)) );
	
	//return 0.5* sqrt( E / (3.0*sumW) );
	//return sqrt( m_E / (3.0*m_sumW) );
}












