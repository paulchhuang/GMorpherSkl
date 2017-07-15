/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#ifndef SOLVER_CERES_H_DEFINED
#define SOLVER_CERES_H_DEFINED

#include "IEnergyTerm_ceres.h"
#include <list>
#include <JointAnglePrior.h>

//#define JOINT_PRIOR

namespace GMorpher 
{


	class SolverCeres
{
	public :
		SolverCeres(const std::vector<GMorpher::rigidT>& RT0,
			const std::vector<CC3D::float3>&	 joint_coords0);
		~SolverCeres();

	// returns the number of iterations performed.. 
	// 0 means we couldnt decrease the energy
	int solve(int                    maxIter,	//ceres::Problem&		problem,
			  const std::list<EPtrCeres>& energies,
	          std::vector<rigidT>&   RT, 
			  std::vector<float3>&	 JX);
	// auxiliary function to transfer the type of JX
	void float3Vector2floatArray(	const std::vector<float3>&	 JX,
									double*   m_JX );
	


	// utility function to compute the energy from outside
	double evaluateEnergy(const std::list<EPtrCeres>& energies,
	                       const std::vector<rigidT>&   RT, 
						   const std::vector<float3>&	JX );

	double evaluateEnergy(const EPtrCeres&              energy,
	                       const std::vector<rigidT>&   RT, 
						   const std::vector<float3>&	JX );
	

	inline int numPatches() const {return m_numPatches; }


#ifdef JOINT_PRIOR
	inline JointAnglePrior*	jointPrior() const{ return m_Jprior; }
#endif

	protected :
	

		int m_numPatches;
		int m_numJoints;
		// solving data
	
		double*             m_RTf64;
		double*				m_JX64;

		unsigned m_concurentThreadsSupported;
		

		// joint angle prior
#ifdef JOINT_PRIOR
	GMorpher::JointAnglePrior*		m_Jprior;
#endif
};




} // end namespace GMorpher













#endif
