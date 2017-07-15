/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#ifndef ENERGYTERM_3DCONSTRAINT_CERES_H_DEFINED
#define ENERGYTERM_3DCONSTRAINT_CERES_H_DEFINED

#include "IEnergyTerm_ceres.h"

namespace GMorpher
{

	struct Point3DConstraint
	{
		Point3DConstraint(const float3* const sourcePoint, const float3* const targetPoint) : m_model(sourcePoint), m_observation(targetPoint)	{}

		template <typename T>
		bool operator()(const T* const para, T* residuals) const {

			T p[3];
			T p0[3];

			p0[0] = T(m_model->x);		p0[1] = T(m_model->y);		p0[2] = T(m_model->z);
			ceres::AngleAxisRotatePoint(para, p0, p);

			p[0] += para[3];
			p[1] += para[4];
			p[2] += para[5];

			residuals[0] = p[0] - T(m_observation->x);
			residuals[1] = p[1] - T(m_observation->y);
			residuals[2] = p[2] - T(m_observation->z);

			return true;
		}

		// Factory to hide the construction of the CostFunction object from the client code.
		static ceres::CostFunction* Create(const float3 const* sourcePoint, const float3 const* targetPoint)
		{
			return (new ceres::AutoDiffCostFunction<Point3DConstraint, 3, 6>(new Point3DConstraint(sourcePoint, targetPoint)));
		}
		const float3* m_model;
		const float3* m_observation;
	};	

class EnergyTerm_3DConstraintCeres : public IEnergyTermCeres
{
	public :
		EnergyTerm_3DConstraintCeres();

		EnergyTerm_3DConstraintCeres(const double  alpha,
							 const double  c,
	                         const int     patch_index,
	                         const float3 const* dX0,
	                         const float3 const* pos);

		virtual ~EnergyTerm_3DConstraintCeres();

		virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX);
		
	inline void updateConstraint( const float3 const* pos) { mPos = pos; }

	protected :
	double          mAlpha;
	double			mC_Tukey;
	int             mPatch_index;
	const float3*   mdX0;
	const float3*   mPos;
};

} 















#endif
