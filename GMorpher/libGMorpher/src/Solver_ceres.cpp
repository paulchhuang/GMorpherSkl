/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#include <Solver_ceres.h>
#include <GraphUtils.h>
#define NUM_JOINT_SUB 16

#include "ceres/ceres.h"
#include "glog/logging.h"

//#define VERBOSE
namespace GMorpher 
{

	SolverCeres::SolverCeres(	const std::vector<GMorpher::rigidT>& RT0,
								const std::vector<CC3D::float3>&	 joint_coords0)
{
	m_numPatches = RT0.size();
	m_numJoints = (joint_coords0.size() > 0) ? joint_coords0.size() - 1 : 0;
	// (jOffsp_bounds.size() > 1) ? jOffsp_bounds.size() - 1 : 0;
	

	// init the solving data
	m_RTf64 = new double[6 * m_numPatches];
	//m_RT_swap = RT0;

	m_JX64 = new double[3 * m_numJoints];
	//m_JX_swap.resize(numJoints);
	m_concurentThreadsSupported = omp_get_num_procs();

#ifdef JOINT_PRIOR	
	m_Jprior	   = new GMorpher::JointAnglePrior("C:\\Users\\Paul\\Desktop\\Ijaz\\pose-conditioned-prior3\\jointAngleModel.mat");
#endif
}


SolverCeres::~SolverCeres()
{
	delete[] m_RTf64;
	delete[] m_JX64;		
	//delete[] m_X;
#ifdef JOINT_PRIOR	
	delete   m_Jprior;
#endif
}





int SolverCeres::solve(int                    maxIter,		//                 ceres::Problem&		problem,
				  const std::list<EPtrCeres>& energies,
                  std::vector<rigidT>&   RT, 
				  std::vector<float3>&	 JX, 
				  bool					 verbose)
{
	
	//assert( RT.size() == m_oadj_bounds.size() - 1);//small consistency check

	typedef std::list<EPtrCeres>::const_iterator EItr;
	int numPatches = RT.size();
	int numJoints = JX.size();

	ceres::Problem		problem;

	// --------------
	// 0 - flatten every unknowns and make them "double"
#pragma omp parallel for
	for (int pi = 0; pi < numPatches; pi++)
	{
		std::vector<float>  angleAxis_tmp(3, .0);
		ceres::QuaternionToAngleAxis<float>((float*)&RT[pi].m_q, angleAxis_tmp.data());

		m_RTf64[6 * pi + 0] = angleAxis_tmp[0];
		m_RTf64[6 * pi + 1] = angleAxis_tmp[1];
		m_RTf64[6 * pi + 2] = angleAxis_tmp[2];
		m_RTf64[6 * pi + 3] = RT[pi].m_t.x;
		m_RTf64[6 * pi + 4] = RT[pi].m_t.y;
		m_RTf64[6 * pi + 5] = RT[pi].m_t.z;
	}
	float3Vector2floatArray(JX, m_JX64);	


	// --------------
	// 1 - add residuals	
	for( EItr e_itr = energies.begin(); e_itr != energies.end(); ++e_itr )
		 (*e_itr)->addResidualBlock(problem, m_RTf64, m_JX64);

	// --------------
	// 2 - solve
	
	ceres::Solver::Options options;
	options.max_num_iterations = maxIter;
	options.num_threads = m_concurentThreadsSupported;		//options.minimizer_type = ceres::LINE_SEARCH;	


	options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.minimizer_progress_to_stdout = false;
	ceres::Solver::Summary summary;

	clock_t start = clock();
	ceres::Solve(options, &problem, &summary);
	clock_t stop = clock();
	
	if (verbose) {
		std::cout << "spend: " << (stop - start) << " ms, " << summary.iterations.size() - 1 << " interations on solving." << std::endl;
		std::cout << summary.FullReport() << "\n";
	}


	// --------------
	// 3 - write the output
#pragma omp parallel for
	for (int pi = 0; pi < numPatches; pi++)
	{
		std::vector<double> q_tmp(4, .0);
		ceres::AngleAxisToQuaternion<double>(m_RTf64 + 6 * pi, q_tmp.data());

		RT[pi].m_q.w = q_tmp[0];
		RT[pi].m_q.x = q_tmp[1];
		RT[pi].m_q.y = q_tmp[2];
		RT[pi].m_q.z = q_tmp[3];
		RT[pi].m_t.x = m_RTf64[6 * pi + 3];
		RT[pi].m_t.y = m_RTf64[6 * pi + 4];
		RT[pi].m_t.z = m_RTf64[6 * pi + 5];
	}


	for (int ji = 0; ji < numJoints; ++ji)
	{
		JX[ji].x = m_JX64[3 * ji + 0];
		JX[ji].y = m_JX64[3 * ji + 1];
		JX[ji].z = m_JX64[3 * ji + 2];		
	}

	

	return (summary.iterations.size() - 1);		// the real number of iterations are iterations.size() - 1

}



void SolverCeres::float3Vector2floatArray(const std::vector<float3>&	 JX,
										double*   m_JX ){
	int numJoint = JX.size();
	double* j_itr = m_JX;
	for(int ji = 0; ji < numJoint; ++ji){
		*(j_itr) = JX[ji].x;		*(j_itr+1) = JX[ji].y;		*(j_itr+2) = JX[ji].z;
		j_itr+=3;
	}
}



} // end namespace GMorpher
