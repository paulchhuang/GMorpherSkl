/* *************************************************
 * Copyright (2012) : Paul Huang, Cedric Cagniart
 * *************************************************/
#ifndef IENERGYTERM_CERES_H_DEFINED
#define IENERGYTERM_CERES_H_DEFINED
#pragma warning(disable: 4819) 
#include "ceres/ceres.h"
#include "ceres/loss_function.h"
#include "ceres/rotation.h"

#include "RigidT.h"
#include <boost\shared_ptr.hpp>


namespace GMorpher
{

	class IEnergyTermCeres 
{
	public :
		virtual ~IEnergyTermCeres();
		virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX) = 0;
	//virtual double computeEnergy( const float* RT, const float* JX) = 0;

	// Functions in the father class cannot overload (best not)
	//virtual void addToGTG_GTb( const float* RT, const float* JX, BigMatrix& BM) const = 0;	// Paul: function overload for bone-binding energy
	//virtual double computeEnergy( const float* RT, const float* JX) const = 0;		// Paul: function overload for bone-binding energy
};

typedef boost::shared_ptr<IEnergyTermCeres> EPtrCeres;


} // end namespace GMorpher 












#endif
