/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#include <EnergyTerm_3DConstraint_ceres.h>

namespace GMorpher
{

	EnergyTerm_3DConstraintCeres::EnergyTerm_3DConstraintCeres() :
		mAlpha(0.),
		mPatch_index(0),
		mC_Tukey(.0)
	{
		
	}

	EnergyTerm_3DConstraintCeres::EnergyTerm_3DConstraintCeres(const double  alpha,
												  const double c,
                                                  const int     patch_index,
                                                  const float3 const* dX0,
                                                  const float3 const* pos):
		mAlpha(alpha),
		mPatch_index(patch_index),
		mdX0(dX0),
		mPos(pos),
		mC_Tukey(c)
	{
	}

	EnergyTerm_3DConstraintCeres::~EnergyTerm_3DConstraintCeres(){}

	void EnergyTerm_3DConstraintCeres::addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX)
	{		
		ceres::CostFunction* cost_function = Point3DConstraint::Create(mdX0, mPos);
		problem.AddResidualBlock(cost_function,
			new ceres::ScaledLoss(new ceres::TukeyLoss(mC_Tukey), mAlpha, ceres::TAKE_OWNERSHIP),/* squared loss */
			RT + 6 * mPatch_index);
	}
	
	//
} // end namespace GMorpher
