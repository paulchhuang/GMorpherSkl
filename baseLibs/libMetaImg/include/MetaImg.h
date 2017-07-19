/* *************************************************
 * Copyright (2017) : Chun-Hao Paul Huang
 * *************************************************/
#ifndef METAIMG_H_DEFINED
#define METAIMG_H_DEFINED

// MVLib
#include <MVLib.h>
// OpenCV
#include <cv.h>
#include <opencv2/highgui/highgui.hpp>


// #############################################################################
// #############################################################################
// MetaImageData
// #############################################################################
// #############################################################################

static void onMouseGRAY(int event, int x, int y, int f, void* image){

	cv::Mat* img = (cv::Mat*)image;

	int w = img->cols;
	int ul = -w - 1; int uc = -w; int ur = -w + 1;
	int cl = -1;                 int cr = -w + 1;
	int ll = w - 1; int lc = w;  int lr = w + 1;

	if (event == cv::EVENT_LBUTTONDOWN){
		unsigned char* pix = (unsigned char*)img->data + y*img->cols + x;

		std::cout << "\n pixel (x, y) = (" << x << "," << y << ") has gray value: " << (int)pix[0] << std::endl;

	}
}


static void onMouseDepth(int event, int x, int y, int f, void* image){

	cv::Mat* img = (cv::Mat*)image;

	int w = img->cols;
	int ul = -w - 1; int uc = -w; int ur = -w + 1;
	int cl = -1;                 int cr = -w + 1;
	int ll = w - 1; int lc = w;  int lr = w + 1;

	if (event == cv::EVENT_LBUTTONDOWN)
	{
		float* pix = (float*)img->data + y*img->cols + x;
		std::cout << "\n pixel (x, y) = (" << x << "," << y << ") has gray value: " << (float)*pix << std::endl;
	}
}


class MetaImgData
{
public:
	
	virtual ~MetaImgData();

	inline  int                  w()      const { return mW; }
	inline  int                  h()      const { return mH; }
	inline  void                 setW(const int w) { mW = w; }
	inline  void                 setH(const int h) { mH = h; }
	inline const double*         GLPMat() const { return mGLPMat; }
	inline const double*         GLMMat() const { return mGLMMat; }
	inline const double*         GLTMat() const { return mGLTMat; }
	inline const double*         P3x4()   const { return mP3x4; }

protected:


	int               mW;
	int               mH;
	double            mGLPMat[16];
	double            mGLMMat[16];
	double            mGLTMat[16];
	double            mP3x4[12];
};



class SilhouetteData : public MetaImgData
{
public:
	SilhouetteData(const MVLib::Camera& cam, const int w, const int h, 	double znear, double zfar);
	virtual ~SilhouetteData();
	
	inline cv::Mat&				 sil()			{ return mSil; }
	

protected:
	cv::Mat			  mSil;		/*[0,255] gray scale image*/

};

#endif
