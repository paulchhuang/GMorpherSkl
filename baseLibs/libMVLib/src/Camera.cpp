/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <MVLib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
namespace MVLib
{


	Matrix3x4 loadPMatrix(const char* filename)
	{
		std::ifstream fin(filename);
		if( !fin.is_open() ) throw std::runtime_error( std::string("could not open :") + std::string(filename) );
		double tempvec[12];
		for(int i=0;i<12;++i) fin >> tempvec[i];
		//Matrix3x4 res(tempvec, tempvec+12);
		Matrix3x4 res;
		res<< tempvec[0],tempvec[1],tempvec[2],tempvec[3],
			  tempvec[4],tempvec[5],tempvec[6],tempvec[7],
			  tempvec[8],tempvec[9],tempvec[10],tempvec[11];
		


		fin.close();
		return res;
	}

	Camera::Camera()
	{
		assert(0 && "we shouldn't be instanciating a camera without matrix");
	}


	Camera::Camera(const Matrix3x4& pmat)
	{
		decompose(pmat, mK, mRw2c, mPos, mPrinRay);
	}

	Camera::Camera(const Camera&B, double xscale, double yscale) :
	mK(B.mK), mPos(B.mPos), mRw2c(B.mRw2c)
	{
		for (int i=0;i<3;++i) {
			mK(0,i) *= xscale;
			mK(1,i) *= yscale;
		}
	}

	// ####################################
	// ACCESSORS
	// ####################################
	const Matrix3& Camera::K() const
	{
		return mK;
	}
	const Matrix3& Camera::Rw2c() const
	{
		return mRw2c;
	}
	const Vector3& Camera::Pos() const
	{
		return mPos;
	}
	const Vector3& Camera::PrinRay() const
	{
		return mPrinRay;
	}
	const Matrix3x4 Camera::P3x4() const
	{
		Matrix3x4 P;
		Matrix3 KR;
		Vector3 KRt;

		/*KR = prod(mK, mRw2c);	
		KRt = prod(KR, mPos);*/		//Paul

		KR = mK*mRw2c;
		KRt = KR*mPos;

		for(int i=0;i<3;++i)
		{
			P(i,3) = - KRt(i);
			for(int j=0;j<3;++j)
				P(i,j) = KR(i,j);
		}
		return P;
	}


	// ####################################
	// EXTRACTORS
	// ####################################
	void Camera::extract_K(double* oK) const
	{
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j)
				oK[3*i+j] = mK(i,j);
		}
	}

	void Camera::extract_Rw2c(double* oRw2c) const
	{
		//std::copy(mRw2c.begin(), mRw2c.end(), oRw2c);
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j)
				oRw2c[3*i+j] = mRw2c(i,j);
		}
	}

	void Camera::extract_Pos(double* oPos) const
	{
		for(int i=0;i<3;++i){
			oPos[i] = mPos(i);
		}
	}

	void Camera::extract_P3x4(double* oP3x4) const
	{
		Matrix3x4 P = P3x4();
		for(int i=0;i<3;++i){
			for(int j=0;j<4;++j){
				oP3x4[4*i+j] = P(i,j);
			}
		}
		//~ for(int i=0;i<3;++i){
			//~ for(int j=0;j<4;++j){
				//~ std::cout<<P(i,j)<<" \t";
			//~ }
			//~ std::cout<<std::endl;
		//~ }
		//~ for(int i=0;i<3;++i){
			//~ for(int j=0;j<4;++j){
				//~ std::cout<<oP3x4[4*i+j]<<" \t";
			//~ }
			//~ std::cout<<std::endl;
		//~ }
	}

	void Camera::extract_GLPROJECTION_Matrix(double* GLPMat, const int width, const int height, const double znear, const double zfar) const
	{
		GLPMat[0*4+0] = (2.0 *  mK(0,0) / double(width));
		GLPMat[1*4+0] = (2.0 *  mK(0,1) / double(width));
		GLPMat[2*4+0] = (((2.0 * mK(0,2) + 1.0)/ double(width))  - 1.0);
		GLPMat[3*4+0] = 0.0;

		GLPMat[0*4+1] = 0.0;
		GLPMat[1*4+1] = (2.0 *  mK(1,1) / double(height));
		GLPMat[2*4+1] = (((2.0 * mK(1,2) + 1.0)/ double(height)) - 1.0);
		GLPMat[3*4+1] = 0.0;

		GLPMat[0*4+2] = 0.0;
		GLPMat[1*4+2] = 0.0;
		GLPMat[2*4+2] = (zfar + znear)/(zfar - znear);
		GLPMat[3*4+2] = -2.0 * zfar * znear / (zfar - znear);

		GLPMat[0*4+3] = 0.0;
		GLPMat[1*4+3] = 0.0;
		GLPMat[2*4+3] = 1.0;
		GLPMat[3*4+3] = 0.0;
	}

	void Camera::extract_GLMODELVIEW_Matrix(double* GLMMat) const
	{
		// 2 -processing
		//Vector3 _Rt( prod( mRw2c,(mPos*-1.0) ) );
		Vector3 ggg = mPos.array()*(-1.0);
		Vector3 _Rt = mRw2c*ggg;

		// 3- converting out
		for (int i=0;i<3;++i)
		{
			for (int j=0;j<3;++j)
			{
				GLMMat[i*4+j] = mRw2c(j,i);
			}
			GLMMat[3*4+i] = _Rt(i);
		}
		GLMMat[0*4+3] = GLMMat[1*4+3] = GLMMat[2*4+3] = 0;
		GLMMat[3*4+3] = 1;
	}

	void Camera::extract_GLTEXTURE2D_Matrix(double* GLT2DMat, const int width, const int height, const double znear, const double zfar) const
	{
		double p[16],m[16],temp[16];
		extract_GLPROJECTION_Matrix(p,width, height, znear, zfar);
		extract_GLMODELVIEW_Matrix(m);

		double scaler[16] = {0.5, 0.0, 0.0, 0.0,
		                     0.0, 0.5, 0.0, 0.0,
		                     0.0, 0.0, 0.5, 0.0,
		                     0.5, 0.5, 0.5, 1.0};

		// HACK : We want to query in the center of pixels. Thus the idea is to add half a pixel to the coordinates
		// since the x and y coordinate will be divided by z, we add 0.5*z to them
		//p[2*4+0] += 1.0/double(width);
		//p[2*4+1] += 1.0/double(height);
							 //~ PRESENTLY THE HACK IS GONE

		//writing the freaking multiplication by hand
		// temp = m * p
		for (int i=0;i<4;++i)
		{
			for (int j=0;j<4;++j)
			{
				temp[i*4+j] = 0;
				for (int k=0;k<4;++k)
				    temp[i*4+j] += m[i*4+k] * p[k*4+j];
			}
		}

		//and another one
		// p = temp * scaler = m * p * scaler
		for (int i=0;i<4;++i)
		{
			for (int j=0;j<4;++j)
			{
				GLT2DMat[i*4+j] = 0;
				for (int k=0;k<4;++k)
					GLT2DMat[i*4+j] += temp[i*4+k] * scaler[k*4+j];
			}
		}
	}








	// ####################################
	// DECOMPOSITION
	// ####################################
	void Camera::decompose(const Matrix3x4& pmat, Matrix3& param, Matrix3& R, Vector3& t, Vector3& prin_ray)
	{
		Vector3 row1, row2, row3, column4;
		//Matrix3 M;
		/*row1    = pmat(0,0), pmat(0,1), pmat(0,2);
		row2    = pmat(1,0), pmat(1,1), pmat(1,2);
		row3    = pmat(2,0), pmat(2,1), pmat(2,2);
		column4 = pmat(0,3), pmat(1,3), pmat(2,3);*/

		row1    << pmat(0,0), pmat(0,1), pmat(0,2);
		row2    << pmat(1,0), pmat(1,1), pmat(1,2);
		row3    << pmat(2,0), pmat(2,1), pmat(2,2);
		column4 << pmat(0,3), pmat(1,3), pmat(2,3);

		/*M  << pmat(0,0), pmat(0,1), pmat(0,2),
		      pmat(1,0), pmat(1,1), pmat(1,2),
		      pmat(2,0), pmat(2,1), pmat(2,2);
*/
		//double gamma = norm2(row3);
		double gamma = row3.norm();

		row1/=gamma;
		row2/=gamma;
		row3/=gamma;
		column4/=gamma;

//		prin_ray = (M.determinant() > 0) ? row3 : -row3;

		//Vector3 q13 = cross(row1,row3);
		//Vector3 q23 = cross(row2,row3);
		Vector3 q13 = row1.cross(row3);
		Vector3 q23 = row2.cross(row3);
		/*double fu = norm2(q13);
		double fv = norm2(q23);*/
		double fu = q13.norm();
		double fv = q23.norm();
		/*double du = dot(row1,row3);
		double dv = dot(row2,row3);*/
		double du = row1.dot(row3);
		double dv = row2.dot(row3);

		/*Vector3 Rrow1( cross(row3,q13) / fu );
		Vector3 Rrow2( cross(row3,q23) / fv );*/
		Vector3 Rrow1( row3.cross(q13) / fu );
		Vector3 Rrow2( row3.cross(q23) / fv );

		R << Rrow1(0),Rrow1(1),Rrow1(2),
			Rrow2(0),Rrow2(1),Rrow2(2),
			row3(0),row3(1),row3(2);

		//param = fu,0,du,0,fv,dv,0,0,1;
		param << fu,0,du,0,fv,dv,0,0,1;


		Vector3 Rt;
		//Rt = (column4(2)*du -column4(0))/fu, (column4(2)*dv -column4(1))/fv , -column4(2);
		Rt << (column4(2)*du -column4(0))/fu, (column4(2)*dv -column4(1))/fv , -column4(2);
		//t = prod(trans(R), Rt);
		t = R.transpose()*Rt;
	}

}
