/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef MVLIB_H_DEFINED
#define MVLIB_H_DEFINED

#include <Eigen/Dense>

#include <list>
#include <map>

namespace MVLib
{

	// #####################################################
	// TYPES
	// #####################################################
	/*typedef tvmet::Matrix<double, 3, 4>	Matrix3x4;
	typedef tvmet::Matrix<double, 3, 3> Matrix3;
	typedef tvmet::Vector<double, 3>    Vector3;*/
	typedef Eigen::Matrix<double, 3, 4>	Matrix3x4;
	typedef Eigen::Matrix<double, 3, 3> Matrix3;
	typedef Eigen::Matrix<double, 3, 1> Vector3;


	class Camera
	{
		public :
		Camera();
		Camera(const Matrix3x4& pmat);
		Camera(const Camera&B, double xscale, double yscale);

		// ---------------------
		// DECOMPOSITION
		//void decompose(const Matrix3x4& pmat, Matrix3& param, Matrix3& R, Vector3& t);
		void decompose(const Matrix3x4& pmat, Matrix3& param, Matrix3& R, Vector3& t, Vector3& prin_ray);
		// ---------------------
		// ACCESSORS
		const Matrix3&   K()    const;
		const Matrix3&   Rw2c() const;
		const Vector3&   Pos()  const;
		const Vector3&   PrinRay()const;
		const Matrix3x4  P3x4() const;
		
		// ---------------------
		// EXTRACTORS
		void extract_K(double* oK)       const;
		void extract_Rw2c(double* oRw2c) const;
		void extract_Pos(double* oPos)   const;
		void extract_P3x4(double* oP3x4) const;
		void extract_PrinRay(double* oPrinRay) const;
		void extract_GLPROJECTION_Matrix(double* GLPMat, const int width, const int height, const double znear, const double zfar) const;
		void extract_GLMODELVIEW_Matrix(double* GLMMat) const;
		void extract_GLTEXTURE2D_Matrix(double* GLT2DMat, const int width, const int height, const double znear, const double zfar) const;

		protected :
		Matrix3 mK;
		Matrix3 mRw2c;
		Vector3 mPos;
		Vector3 mPrinRay;
	};

	// #####################################################
	// FUNCTIONS
	// #####################################################
	std::string buildFilename(const char* baseName, int id);
	Matrix3x4 loadPMatrix(const char* filename);


	// ----------------------------------
	// CameraSetup.cpp
	void loadCameras(const char*            cameraCalibBaseName,
	                 const std::list<int>&  indexList,
	                 std::map<int, Camera>& cams, 
					 const bool KinovisFlag = false);

	void getNearFar(const std::map<int, Camera>& cams,
	                const double xmin, const double ymin, const double zmin,
	                const double xmax, const double ymax, const double zmax,
	                double& znear, double& zfar);

	
	// ----------------------------------
	// ColorConvert.cpp
	void RGB_2_GRAY( const unsigned char* RGB_buffer, unsigned char* GRAY_buffer, const int w, const int h);
	void GRAY_2_RGB( const unsigned char* GRAY_buffer, unsigned char* RGB_buffer, const int w, const int h);

	// ----------------------------------
	// SubSampler.cpp
	void subsampleN_gray8(const int N,
	                      unsigned char* bmp_src, int w_src, int h_src,
						  unsigned char* bmp_dst, int w_dst, int h_dst,
						  unsigned char* WS = 0);
	void subsampleN_rgb8(const int N,
	                      unsigned char* bmp_src, int w_src, int h_src,
						  unsigned char* bmp_dst, int w_dst, int h_dst,
						  unsigned char* WS = 0);
	void subsampleN_rgba8(const int N,
	                      unsigned char* bmp_src, int w_src, int h_src,
	                      unsigned char* bmp_dst, int w_dst, int h_dst,
	                      unsigned char* WS = 0);
	void subsampleN_gray32(const int N,
	                      unsigned char* bmp_src, int w_src, int h_src,
						  unsigned char* bmp_dst, int w_dst, int h_dst,
						  unsigned char* WS = 0);

	// ----------------------------------
	// LineDrawer.cpp
	void drawLineGRAY( unsigned char* image, int W, int x1, int y1, int x2, int y2, unsigned char val );
	void drawLineRGB( unsigned char* image, int W, int x1, int y1, int x2, int y2, unsigned char r, unsigned char g, unsigned char b );

	void findEdges(const int w, const int h, unsigned char* const src, unsigned char* dst);
	void findEdges(const int w, const int h, unsigned char* const src, bool* dst);

	/**
	* Auxilliary functor...
	* takes an image as input, and computes the gradient where it is requested
	*/
	template <typename T>
	class GradientCompute
	{
		public :
		inline GradientCompute(T* image, int w)
		{
			_image = image;
			dmm = -w-1; dcm = -w; dpm = -w +1;
			dmc = -1;             dpc = 1;
			dmp = w-1;  dcp = w;  dpp = w +1;
		}

		inline void operator () (int px_offset, double& gx, double& gy)
		{// sobel mask
			T* ptr = _image + px_offset;

			gx =-(0.5*(ptr[dmm] - ptr[dpm]));
			gx -=     (ptr[dmc] - ptr[dpc]);
			gx -= 0.5*(ptr[dmp] - ptr[dpp]);

			gy =-( ptr[dcm] + 0.5*(ptr[dmm] + ptr[dpm]));
			gy +=  ptr[dcp] + 0.5*(ptr[dmp] + ptr[dpp]);
		}

		protected :
		T*  _image;
		int dmm, dcm, dpm;
		int dmc,      dpc;
		int dmp, dcp, dpp;
	};


	template <typename T>
	class RGBGradientCompute
	{
		public :
		inline RGBGradientCompute(T* image, int w)
		{
			_image = image;
			dmm = 3*(-w-1); dcm = -3*w; dpm = 3*(-w +1);
			dmc = -3;					dpc = 3;
			dmp = 3*(w-1);  dcp = 3*w;	dpp = 3*(w +1);
		}

		inline double operator () (int px_offset, double& gx, double& gy)
		{// sobel mask
			T* ptr = _image + px_offset;		// note: here the offset is already 3*(row*w + col)

			double gx_tmp, gy_tmp, magSQ_tmp, MAX_magSQ = 0;

			for (int c = 0; c < 3; ++c){
				gx_tmp =-(0.5*(ptr[dmm] - ptr[dpm]));
				gx_tmp -=     (ptr[dmc] - ptr[dpc]);
				gx_tmp -= 0.5*(ptr[dmp] - ptr[dpp]);

				gy_tmp =-( ptr[dcm] + 0.5*(ptr[dmm] + ptr[dpm]));
				gy_tmp +=  ptr[dcp] + 0.5*(ptr[dmp] + ptr[dpp]);

				magSQ_tmp = gx_tmp*gx_tmp + gy_tmp*gy_tmp;

				if (magSQ_tmp > MAX_magSQ){
					gx			= gx_tmp;		
					gy			= gy_tmp;
					MAX_magSQ	= magSQ_tmp;
				}
				ptr++;
			}

			return sqrt(MAX_magSQ);
		}

		protected :
		T*  _image;
		int dmm, dcm, dpm;
		int dmc,      dpc;
		int dmp, dcp, dpp;
	};

}


#endif