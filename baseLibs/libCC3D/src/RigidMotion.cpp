/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <CC3D/RigidMotion.h>
#include <iostream>
#include <cmath>

extern "C" {
	void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
	void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
		                  double *A, int *LDAp, double *VLp, double *VUp,
		                  int *ILp, int *IUp, double *ABSTOLp, int *Mp,
		                  double *W, double *Z, int *LDZp, int *ISUPPZ,
		                  double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
		                  int *INFOp);

}


namespace CC3D
{

	// -------------------------------------------------
	// -------------------------------------------------
	// THE POLAR DECOMPOSITION
	// -------------------------------------------------
	// -------------------------------------------------
	int RigidMotionPolar::LWORK = 18;
	char RigidMotionPolar::JOBU = 'A'; //all M columns of U are returned in array U
	char RigidMotionPolar::JOBVT = 'A'; // all N rows of V**T are returned in the array VT;  (we need the last column)
	int RigidMotionPolar::M = 3;
	int RigidMotionPolar::N = 3;
	int RigidMotionPolar::LDU = 3;
	int RigidMotionPolar::LDA = 3;
	int RigidMotionPolar::LDVT = 3;

	int RigidMotionPolar::operator () (const double* covMat, double* R)
	{
		double covMatLocal[9];
		std::copy(covMat, covMat+9, covMatLocal);

		INFO = 0;
		dgesvd_(&JOBU, &JOBVT,
				&M, &N,
				covMatLocal, &LDA,
				S,
				U, &LDU,
				VT, &LDVT,
				Work,&LWORK,
				&INFO);


		// -----------------------------------------------------------------------------------------------------------------------------
		// 3 - write the result to the Rotation Matrix
		// we fed cov = E[x1 x2t] in fortran who thought it was E[x2 x1t]
		// we get cov = (Vft)'*Sf*Uf
		// so if cov = U S Vt
		// we have U = Vft and V = Uft
		// so if R = V* Ut
		// R = Uft * Vf
		// the swapping of the last column of U becomes the swapping of the last column of Vft

		for(int r=0; r<3; ++r) {
			for(int c=0; c<3; ++c) {
				double accum = 0.0;
				for(int k=0;k<3;++k) {
					accum += U[3*k+r] * VT[3*c+k];
				}
				R[3*r+c] = accum;
			}
		}

		double det =   R[3*0+0]* ( R[3*1+1] * R[3*2+2] - R[3*2+1] * R[3*1+2] )
					 - R[3*1+0]* ( R[3*0+1] * R[3*2+2] - R[3*2+1] * R[3*0+2] )
					 + R[3*2+0]* ( R[3*0+1] * R[3*1+2] - R[3*1+1] * R[3*0+2] );

		// flip the last column of U (the smallest eigenvalue) if we have a non special orthogonal matrix
		if(det < 0.0 ) {
			VT[2] = -VT[2]; VT[5] = - VT[5]; VT[8] = - VT[8];
			for(int r=0; r<3; ++r) {
				for(int c=0; c<3; ++c) {
					double accum = 0.0;
					for(int k=0;k<3;++k) {
						accum += U[3*k+r] * VT[3*c+k];
					}
					R[3*r+c] = accum;
				}
			}
		}

		return INFO;
	}
















	// -------------------------------------------------
	// -------------------------------------------------
	// THE HORN DECOMPOSITION
	// -------------------------------------------------
	// -------------------------------------------------
	char    RigidMotionHorn::JOBZ = 'V'; //Compute eigenvalues and eigenvectors.
	char    RigidMotionHorn::RANGE = 'I'; //the IL-th through IU-th eigenvalues will be found.
	char    RigidMotionHorn::UPLO = 'U'; //  Upper triangle of A is stored;
	int     RigidMotionHorn::N = 4;
	int     RigidMotionHorn::LDA = 4;
	double  RigidMotionHorn::VL = 0; //we dont need these (used for 'V' option )
	double  RigidMotionHorn::VU = 0; //we dont need these (used for 'V' option )
	int     RigidMotionHorn::IL = 4;
	int     RigidMotionHorn::IU = 4;
	double  RigidMotionHorn::ABSTOL = 1e-4; // tolerance
	int     RigidMotionHorn::LDZ = 4;
	int     RigidMotionHorn::LWORK = 26*4; //HACK
	int     RigidMotionHorn::LIWORK = 10*4;

	int RigidMotionHorn::operator () (double* covMat, double* R)
	{
		double qw,qx,qy,qz;
		int INFO = getQuaternion( covMat, qw,qx,qy,qz );
		R[0] = 1 - 2*(qy*qy + qz*qz) ;
		R[1] = 2* (qx*qy - qz*qw);
		R[2] = 2* (qx*qz + qy*qw);

		R[3] = 2* (qx*qy + qz*qw);
		R[4] = 1 - 2*(qx*qx + qz*qz) ;
		R[5] = 2* (qy*qz - qx*qw);

		R[6] = 2* (qx*qz - qy*qw);
		R[7] = 2* (qy*qz + qx*qw);
		R[8] = 1 - 2*(qx*qx + qy*qy) ;

		return INFO;
	}



	int RigidMotionHorn::getQuaternion( const double* covMat, double& qw, double& qx, double& qy, double& qz )
	{
		double Nmat[16];
		// diagonal terms first
		Nmat[4*0 + 0] =   covMat[0*3+0] + covMat[1*3+1] + covMat[2*3+2]; //  cxx + cyy +czz
		Nmat[4*1 + 1] =   covMat[0*3+0] - covMat[1*3+1] - covMat[2*3+2]; //  cxx - cyy -czz
		Nmat[4*2 + 2] = - covMat[0*3+0] + covMat[1*3+1] - covMat[2*3+2]; // -cxx + cyy -czz
		Nmat[4*3 + 3] = - covMat[0*3+0] - covMat[1*3+1] + covMat[2*3+2]; // -cxx - cyy + czz
		// off diagonal
		Nmat[4*0 + 1] = Nmat[0 + 4*1] = covMat[3*1+2] - covMat[3*2+1]; //cyz - czy
		Nmat[4*0 + 2] = Nmat[0 + 4*2] = covMat[3*2+0] - covMat[3*0+2]; //czx - cxz
		Nmat[4*0 + 3] = Nmat[0 + 4*3] = covMat[3*0+1] - covMat[3*1+0]; //cxy - cyx
		Nmat[4*1 + 2] = Nmat[1 + 4*2] = covMat[3*0+1] + covMat[3*1+0]; //cxy + cyx
		Nmat[4*1 + 3] = Nmat[1 + 4*3] = covMat[3*2+0] + covMat[3*0+2]; //czx + cxz
		Nmat[4*2 + 3] = Nmat[2 + 4*3] = covMat[3*1+2] + covMat[3*2+1]; //cyz + czy


		int    M = 0 ;//number of eigen values found
		double W[4]; // storage for eigenvalues
		double Z[16]; // storage for eigenvectors
		int    ISUPPZ[8];
		double WORK[26*4]; //HACK

		int IWORK[10*4];  //HACK
		int INFO = 0;

		dsyevr_( &JOBZ, &RANGE, &UPLO, &N,
		         Nmat, &LDA, &VL, &VU,
		         &IL, &IU, &ABSTOL, &M,
		         W, Z, &LDZ, ISUPPZ,
		         WORK, &LWORK, IWORK, &LIWORK,
		         &INFO);

		qw = Z[0];
		qx = Z[1];
		qy = Z[2];
		qz = Z[3];

		return INFO;
	}
























	// -------------------------------------------------
	// -------------------------------------------------
	// THE POLAR DECOMPOSITION ( MY VERSION )
	// -------------------------------------------------
	// -------------------------------------------------

	/**
	 * returns such that A = Q diag QT
	 */
	int myjacobi(const double* A, double* diag, double* Q)
	{


		double offd[3];

		Q[3*0+0]=1.; Q[3*0+1]=0.; Q[3*0+2]=0.;
		Q[3*1+0]=0.; Q[3*1+1]=1.; Q[3*1+2]=0.;
		Q[3*2+0]=0.; Q[3*2+1]=0.; Q[3*2+2]=1.;

		diag[0] = A[0*3+0]; diag[1] = A[1*3+1]; diag[2] = A[2*3+2]; // copy the diag part
		offd[0] = A[1*3+2]; offd[1] = A[2*3+0]; offd[2] = A[0*3+1]; // copy the offdiag part


		int sweep = 0 ;
		for(;sweep<20; ++sweep)
		{
			double sum_off_diag = fabs(offd[0]) + fabs(offd[1]) + fabs(offd[2]);
			if (sum_off_diag == 0.0f) break;

			for( int i=0;i<3;++i) {
				 // first we find the two other coordinates on which we ll operate
				int p = (i+1)%3;
				int q = (p+1)%3;
				double fabsOffDi = fabs( offd[i] ); // this is the a[p][q]
				double g = 100.0f*fabsOffDi;
				if( fabsOffDi > 0.0f) {
					double h = diag[q] - diag[p];
					double fabsh = fabs(h);
					double t;
					if( fabsh+g == fabsh){
						t = offd[i]/h;
					}
					else {
						double theta = 0.5f*h/offd[i];
						t = 1.0f/(fabs(theta) + sqrt(1.0f+theta*theta));
						if ( theta < 0. ) t = -t;
					}
					double c = 1.0f/sqrt(1.0f+t*t);
					double s   = t*c;
					double tau = s/(1.0f+c);
					double ta   = t*offd[i];
					offd[i] = 0.0f;

					// updating the eigen values
					diag[p] -= ta;
					diag[q] += ta;

					// updating the matrix
					double offdq = offd[q]; // a swap term
					offd[q] -= s*(offd[p] + tau*offd[q]);
					offd[p] += s*(offdq   - tau*offd[p]);
					for( int j=0;j<3;++j) {
						double Qjp = Q[3*j+p];
						double Qjq = Q[3*j+q];
						Q[3*j+p] -= s*(Qjq + tau*Qjp);
						Q[3*j+q] += s*(Qjp - tau*Qjq);
					}
				}
			}
		}


		int I1 = 0;
		if( diag[1] > diag[0] )  I1 = 1;
		if( diag[2] > diag[I1] ) I1 = 2;
		int I2 = (I1+1)%3;
		int I3 = (I1+2)%3;
		if( diag[I3] > diag[I2] ) {int t=I2; I2 = I3; I3 = t; }

		if( diag[0] < 0 ) diag[0] = 0.;
		if( diag[1] < 0 ) diag[1] = 0.;
		if( diag[2] < 0 ) diag[2] = 0.;


		double Qswap[9];
		for(int i=0;i<9;++i) Qswap[i] = Q[i];
		Q[0] = Qswap[I1]; Q[3] = Qswap[I1+3]; Q[6] = Qswap[I1+6];
		Q[1] = Qswap[I2]; Q[4] = Qswap[I2+3]; Q[7] = Qswap[I2+6];
		Q[2] = Qswap[I3]; Q[5] = Qswap[I3+3]; Q[8] = Qswap[I3+6];

		Qswap[0] = diag[I1];
		Qswap[1] = diag[I2];
		Qswap[2] = diag[I3];

		diag[0] = Qswap[0];
		diag[1] = Qswap[1];
		diag[2] = Qswap[2];
		return sweep;
	}





	int RigidMotionPolar2::operator () (const double* covMat, double* R)
	{
		// build ATA
		double ATA[9];
		for(int i=0;i<3;++i) {
			for(int j=0;j<3;++j) {
				ATA[3*i+j] = covMat[3*0+i]*covMat[3*0+j] + covMat[3*1+i]*covMat[3*1+j] + covMat[3*2+i]*covMat[3*2+j];
			}
		}

		double U[9], S2[3], V[9];
		myjacobi(ATA, S2, V);

		// if we dont even have a proper sing value the mat is 0 and we should give up
		if (S2[0] < 1e-10) return -1;
		// count rank.. it it's 1 we ll have bullshit ( ambiguity in one plane rotation) --> give up
		int rank = 3;
		for( int i=0; i<3; ++i ) if( S2[i]/S2[0] < 1e-6 ) rank--;
		if( rank < 2 ) return -1;

		// copy the full rank part
		for(int j=0; j<rank;++j) {
			for( int i=0;i<3; ++i) U[3*i+j] = ( covMat[3*i+0]*V[3*0+j]
											  + covMat[3*i+1]*V[3*1+j]
											  + covMat[3*i+2]*V[3*2+j] ) / sqrt(S2[j]);
		}

		// for the kernel of AT we need to build the last guy using the first two
		if( rank == 2 ) {
			U[3*0+2] = U[3*1+0] * U[3*2+1] - U[3*2+0] * U[3*1+1];
			U[3*1+2] = U[3*2+0] * U[3*0+1] - U[3*0+0] * U[3*2+1];
			U[3*2+2] = U[3*0+0] * U[3*1+1] - U[3*1+0] * U[3*0+1];
		}


		// ------------------
		// 3 - R = V UT
		for(int i=0;i<3;++i) {
			for(int j=0;j<3;++j) {
				R[3*i+j] = V[3*i+0] * U[3*j+0]
						 + V[3*i+1] * U[3*j+1]
						 + V[3*i+2] * U[3*j+2];
			}
		}


		double det =   R[3*0+0]* ( R[3*1+1] * R[3*2+2] - R[3*2+1] * R[3*1+2] )
					 - R[3*1+0]* ( R[3*0+1] * R[3*2+2] - R[3*2+1] * R[3*0+2] )
					 + R[3*2+0]* ( R[3*0+1] * R[3*1+2] - R[3*1+1] * R[3*0+2] );

		// if non SO3 flip the last column of U (smallest singular val)
		if(det < 0.0 ) {
			U[2] = -U[2]; U[5] = -U[5]; U[8] = -U[8];
			for(int i=0;i<3;++i) {
				for(int j=0;j<3;++j) {
				R[3*i+j] = V[3*i+0] * U[3*j+0]
						 + V[3*i+1] * U[3*j+1]
						 + V[3*i+2] * U[3*j+2];
				}
			}
		}


		return 0;
	}

}
