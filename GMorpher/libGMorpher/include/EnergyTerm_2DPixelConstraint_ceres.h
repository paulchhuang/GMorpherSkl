/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#ifndef ENERGYTERM_2DPIXELCONSTRAINT_CERES_H_DEFINED
#define ENERGYTERM_2DPIXELCONSTRAINT_CERES_H_DEFINED

#include "IEnergyTerm_ceres.h"

namespace GMorpher
{

	struct Pixel2DConstraint
	{
		Pixel2DConstraint(const float3* const sourcePoint, const std::pair<int, int>& targetPixel, const double* const P3x4) : m_model(sourcePoint), m_utarget(targetPixel.first), m_vtarget(targetPixel.second), m_P3x4(P3x4){}

		template <typename T>
		bool operator()(const T* const para, T* residuals) const {

			T p[3];
			T p0[3];

			p0[0] = T(m_model->x);		p0[1] = T(m_model->y);		p0[2] = T(m_model->z);
			ceres::AngleAxisRotatePoint(para, p0, p);

			p[0] += para[3];
			p[1] += para[4];
			p[2] += para[5];

			T xcam = m_P3x4[0 * 4 + 0] * p[0] + m_P3x4[0 * 4 + 1] * p[1] + m_P3x4[0 * 4 + 2] * p[2] + m_P3x4[0 * 4 + 3];
			T ycam = m_P3x4[1 * 4 + 0] * p[0] + m_P3x4[1 * 4 + 1] * p[1] + m_P3x4[1 * 4 + 2] * p[2] + m_P3x4[1 * 4 + 3];
			T zcam = m_P3x4[2 * 4 + 0] * p[0] + m_P3x4[2 * 4 + 1] * p[1] + m_P3x4[2 * 4 + 2] * p[2] + m_P3x4[2 * 4 + 3];

			T ucam = xcam / zcam;
			T vcam = ycam / zcam;

			residuals[0] = ucam - T(m_utarget);
			residuals[1] = vcam - T(m_vtarget);

		
			return true;
		}

		// Factory to hide the construction of the CostFunction object from the client code.
		static ceres::CostFunction* Create(const float3 const* sourcePoint, const std::pair<int, int>& targetPixel, const double* const P3x4)
		{
			return (new ceres::AutoDiffCostFunction<Pixel2DConstraint, 2, 6>(new Pixel2DConstraint(sourcePoint, targetPixel, P3x4)));
		}
		const float3* m_model;
		const double* m_P3x4;
		const int m_utarget;
		const int m_vtarget;
	};	

//class EnergyTerm_3DConstraintCeres : public IEnergyTermCeres
//{
//	public :
//		EnergyTerm_3DConstraintCeres();
//
//		EnergyTerm_3DConstraintCeres(const double  alpha,
//							 const double  c,
//	                         const int     patch_index,
//	                         const float3 const* dX0,
//	                         const float3 const* pos);
//
//		virtual ~EnergyTerm_3DConstraintCeres();
//
//		virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX);
//		
//	inline void updateConstraint( const float3 const* pos) { mPos = pos; }
//
//	protected :
//	double          mAlpha;
//	double			mC_Tukey;
//	int             mPatch_index;
//	const float3*   mdX0;
//	const float3*   mPos;
//};

} 















#endif
