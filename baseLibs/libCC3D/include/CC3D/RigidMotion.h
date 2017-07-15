/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef CC3D_RIGIDMOTION_H_DEFINED
#define CC3D_RIGIDMOTION_H_DEFINED

#include <algorithm>

namespace CC3D
{

	// -------------------------------------------------
	// -------------------------------------------------
	// THE COVMATBUILDER CLASS : builds a covariance matrix incrementally
	// -------------------------------------------------
	// -------------------------------------------------
	class CovMatBuilder
	{
		public :
		inline CovMatBuilder();
		inline void pushCorrespondance(double x1, double y1, double z1,
		                               double x2, double y2, double z2,
		                               double w);
		inline int getCovMat( double* o_covMat, double* o_mean1, double* o_mean2);

		protected :
		double covMat[9];
		double mean1[3];
		double mean2[3];
		double sumW;
	};


	CovMatBuilder::CovMatBuilder()
	{
		std::fill(covMat,covMat+9,0);
		std::fill(mean1, mean1+3,0);
		std::fill(mean2, mean2+3,0);
		sumW = 0;
	}

	void CovMatBuilder:: pushCorrespondance(double x1, double y1, double z1,
		                                    double x2, double y2, double z2,
		                                    double w)
	{
		covMat[3*0+0] += x1 * x2 * w;	covMat[3*0+1] += x1 * y2 * w;	covMat[3*0+2] += x1 * z2 * w;
		covMat[3*1+0] += y1 * x2 * w;	covMat[3*1+1] += y1 * y2 * w;	covMat[3*1+2] += y1 * z2 * w;
		covMat[3*2+0] += z1 * x2 * w;	covMat[3*2+1] += z1 * y2 * w;	covMat[3*2+2] += z1 * z2 * w;

		mean1[0] += x1 * w;		mean1[1] += y1 * w;		mean1[2] += z1 * w;
		mean2[0] += x2 * w;		mean2[1] += y2 * w;		mean2[2] += z2 * w;

		sumW += w;
	}

	int CovMatBuilder::getCovMat( double* o_covMat, double* o_mean1, double* o_mean2)
	{
		if(sumW == 0.0 ) return -1;

		// recover the translation
		for(int i=0; i<3; ++i) {
			o_mean1[i] = mean1[i] /sumW;
			o_mean2[i] = mean2[i] /sumW;
		}

		// recover the covariance matrix
		for(int i=0; i<9 ; ++i)
			o_covMat[i] = covMat[i] / sumW;

		for(int r = 0; r < 3 ; ++r) {
			for(int c = 0; c < 3; ++c) o_covMat[3*r+c] -= o_mean1[r] * o_mean2[c];
		}

		return 0;
	}

	// -------------------------------------------------
	// -------------------------------------------------
	// THE POLAR DECOMPOSITION
	// -------------------------------------------------
	// -------------------------------------------------
	class RigidMotionPolar
	{
		public :

		int operator () (const double* covMat, double* R);


		private :

		double U[9];
		double VT[9];
		double S[3];
		double Work[18];	//the size should be 4*min(M,N) + 2*max(M,N), so it is 6*3 in our case
		int INFO;

		static  int LWORK;
		static  char JOBU;
		static  char JOBVT;
		static  int M;
		static  int N;
		static  int LDU;
		static  int LDA;
		static  int LDVT;
	};

	// -------------------------------------------------
	// -------------------------------------------------
	// THE POLAR DECOMPOSITION ( MY VERSION )
	// -------------------------------------------------
	// -------------------------------------------------
	class RigidMotionPolar2 {
		public :
		int operator () (const double* covMat, double* R);
	};


	// -------------------------------------------------
	// -------------------------------------------------
	// THE HORN DECOMPOSITION
	// -------------------------------------------------
	// -------------------------------------------------
	class RigidMotionHorn
	{
		public :

		int operator () (double* covMat, double* R);
		int getQuaternion(const double* covMat, double& qw, double& qx, double& qy, double& qz );

		private :
		static char JOBZ, RANGE, UPLO;
		static int  N, LDA, IL, IU, LDZ, LWORK, LIWORK;
		static double ABSTOL, VL, VU;
	};


}






#endif
