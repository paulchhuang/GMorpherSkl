/* *************************************************
 * Copyright (2017) : Chun-Hao Paul Huang
 * *************************************************/

#include<MetaImg.h>

using namespace cv;

SilhouetteData::SilhouetteData(const MVLib::Camera& cam,
                                             const int w, const int h,
                                             double znear, double zfar)
{
	this->setH(h);
	this->setW(w);
	cam.extract_GLPROJECTION_Matrix( mGLPMat, w, h, znear, zfar);
	cam.extract_GLMODELVIEW_Matrix ( mGLMMat );
	cam.extract_GLTEXTURE2D_Matrix ( mGLTMat, w, h, znear, zfar);
	cam.extract_P3x4( mP3x4);	
}

MetaImgData::~MetaImgData()
{
	//delete [] mSil;
}

SilhouetteData::~SilhouetteData()
{
	//delete [] mSil;
}
