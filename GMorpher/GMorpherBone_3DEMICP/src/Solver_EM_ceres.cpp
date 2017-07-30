/* *************************************************
 * Copyright (2016) : Paul Huang
 * *************************************************/
#include <Solver_EM_ceres.h>
#include <GraphUtils.h>
#include <numeric>

//#define VERBOSE
Solver_EMCeres::Solver_EMCeres(bool                                        prediction,
								const SklRgedPatchedCloud::PatchedCloud&    PC)

{
	m_GSolver = new GMorpher::SolverCeres(PC.RT0(), PC.joint_coords0());// , sadj, sadj_bounds, jOffsp, jOffsp_bounds);

	if( prediction ) m_EICP = boost::shared_ptr<IEnergyTerm_EMCeres>( 
	                                          new EnergyTerm_EMPredCeres(10, PC) );
	/*else             m_EICP = boost::shared_ptr<IEnergyTerm_EMCeres>( 
	                                              new EnergyTerm_EMCeres(10, PC) );	*/
}


Solver_EMCeres::~Solver_EMCeres()
{
}


void Solver_EMCeres::setRefMesh(SklRgedPatchedCloud::PatchedCloud&  mesh){
	m_EICP->reinit(mesh);
}


bool countMoversCeres( const std::vector<GMorpher::rigidT>& RT_old,
               const std::vector<GMorpher::rigidT>& RT_new,
               const double meanEdge )
{
	double thresh = 0.005*meanEdge*meanEdge;
	// we can evaluate how many centers have moved
	int numPatches = RT_old.size();
	int numMovers = 0;

	std::vector<int> movers(numPatches);

	#pragma omp parallel for 
	for(int pi=0;pi<numPatches;++pi) {
		GMorpher::float3 delta = RT_new[pi].m_t - RT_old[pi].m_t;
		if( dot(delta,delta) > thresh ) movers[pi]++;
		//if( dot(delta,delta) > thresh ) numMovers++;		
	}
	// we can also evaluate how many center predictions have moved ?
	// TODO 	
	return std::accumulate(movers.begin(),movers.end(),0);	
} 


int Solver_EMCeres::solve(const std::list<GMorpher::EPtrCeres>& regularizationTerms,					 
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
					  std::vector<CC3D::float3>&	   JX,	
					  float&                           sigma,
					  bool							   verbose)
{
	int numPatches = m_GSolver->numPatches();

	float* RTf = new float[12*numPatches];
	std::list<GMorpher::EPtrCeres> energies = regularizationTerms;
	energies.push_back(m_EICP);

	std::vector<GMorpher::rigidT> RT_old = RT;
	std::vector<CC3D::float3>	  JX_old = JX;
	double sigmaStart = sigma;
	// set the data
	m_EICP->setFrameData( obs_cloud_coords, obs_cloud_normals);
	m_EICP->setAlpha( alphaNN );

	int   numIter_EM = 0;
	RigidTFlatten( RT, RTf );

	
	bool stop_flag = false;
	std::vector<int>			isValid_patch(numPatches, 1);


	while ( numIter_EM++ < maxIter_EM)
	{
		// E - Step
		if( !probabilistic ) 
		{
			m_EICP->EStep( sigmaStart, normThresh, EOutlier, RTf );			
			m_EICP->EStep_makeDeterministicAssignment();	// if we are using a non-probabilistic approach
		}
		else 
		{			
			m_EICP->EStep( sigma, normThresh, EOutlier, RTf );		
		}

		// M - Step
		if (verbose) std::cout << "# EM_iter.: " << numIter_EM << " ";
		m_GSolver->solve(maxIter_M, energies, RT, JX, verbose);

		// re-evaluate sigma
		RigidTFlatten( RT, RTf );
		sigma = m_EICP->reEvaluateSigma( RTf );
		// check for convergence
		int numMoved = countMoversCeres(RT_old, RT, meanEdge);

		#ifdef VERBOSE 
			std::cout<<"sigma "<<sigma/meanEdge<<" and movers "<<numMoved<<std::endl;
		#endif;
		
		if( numMoved == 0 ) break;
		RT_old = RT;
		JX_old = JX;

		
	}

	delete[] RTf;
	return ((numIter_EM > maxIter_EM) ? (numIter_EM-1) : numIter_EM);
}





