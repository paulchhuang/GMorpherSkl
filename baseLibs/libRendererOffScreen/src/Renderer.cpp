/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
//#define DEBUG_OPENGL
#include <Renderer.h>
#include "GLerrors.h"
#include <iostream>
#include <cassert>

namespace OffScreenRendering
{

// ###################################################
// ###################################################
// RENDERER2D
// ###################################################
// ###################################################

Renderer::Renderer(unsigned int maxW, unsigned int maxH,
				   ColorFormats cfmt, DepthFormats dfmt)
{
	glGenTextures(1, &mColorTex);
	glGenTextures(1, &mDepthTex);
	glGenFramebuffers(1,&mGLFBOId);
	resize(maxW, maxH, cfmt, dfmt);
}

Renderer::~Renderer()
{
	//delete FBO
	glDeleteFramebuffers(1,&mGLFBOId);
	glDeleteTextures(1, &mColorTex);
	glDeleteTextures(1, &mDepthTex);

	checkGLErrors("OFFSCREENRENDERER : destructor");
}



void Renderer::draw(IGLRenderJob* job,
					unsigned int x, unsigned int y, unsigned int dx, unsigned int dy,
					unsigned char* out_image,
					float*         out_depth)
{
//	std::cout<<"rendering "<<dx<<"x"<<dy<<std::endl;
	assert( (x+dx <= mMaxWidth) && (y+dy <= mMaxHeight) );

	// bind the render target
	glBindFramebuffer(GL_FRAMEBUFFER, mGLFBOId);
	// run the rendering job
	job->run();

	// copy the data
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(x, y, dx, dy, mColor_format, mColor_type, out_image);
	glReadPixels(x, y, dx, dy, mDepth_format, mDepth_type, out_depth);

	// return to normal rendering
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	checkGLErrors("OFFSCREENRENDERER : end draw");
}



void Renderer::draw(IGLRenderJob* job,
					unsigned int x, unsigned int y, unsigned int dx, unsigned int dy,
					float*         out_image,
					float*         out_depth)
{
	std::cout<<"rendering "<<dx<<"x"<<dy<<std::endl;
	assert( (x+dx <= mMaxWidth) && (y+dy <= mMaxHeight) );

	// bind the render target
	glBindFramebuffer(GL_FRAMEBUFFER, mGLFBOId);
	// run the rendering job
	job->run();

	// copy the data
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(x, y, dx, dy, mColor_format, mColor_type, out_image);
	glReadPixels(x, y, dx, dy, mDepth_format, mDepth_type, out_depth);

	// return to normal rendering
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	checkGLErrors("OFFSCREENRENDERER : end draw");
}


void Renderer::draw(IGLRenderJob* job,
					unsigned int x, unsigned int y, unsigned int dx, unsigned int dy,
					unsigned char* out_image,
					float*         out_depth,
					unsigned char* out_stencil)
{
	std::cout<<"rendering "<<dx<<"x"<<dy<<std::endl;
	assert( (x+dx <= mMaxWidth) && (y+dy <= mMaxHeight) );

	// bind the render target
	glBindFramebuffer(GL_FRAMEBUFFER, mGLFBOId);
	// run the rendering job
	job->run();

	// copy the data
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(x, y, dx, dy, mColor_format, mColor_type, out_image);
	glReadPixels(x, y, dx, dy, GL_DEPTH_COMPONENT, GL_FLOAT, out_depth);
	glReadPixels(x, y, dx, dy, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, out_stencil);

	// return to normal rendering
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	checkGLErrors("OFFSCREENRENDERER : end draw");
}



void Renderer::resize(unsigned int maxW, unsigned int maxH,
					  ColorFormats cfmt, DepthFormats dfmt)
{
	switch(cfmt)
	{
		case RGB8  : mColor_internalFormat = GL_RGB8;              mColor_format = GL_RGB;       mColor_type = GL_UNSIGNED_BYTE; break;
		case RGBA8 : mColor_internalFormat = GL_RGBA8;             mColor_format = GL_RGBA;      mColor_type = GL_UNSIGNED_BYTE; break;
		case GRAY8 : mColor_internalFormat = GL_LUMINANCE8;        mColor_format = GL_LUMINANCE; mColor_type = GL_UNSIGNED_BYTE; break;
		case GRAY32 : mColor_internalFormat = GL_LUMINANCE32F_ARB; mColor_format = GL_LUMINANCE; mColor_type = GL_FLOAT;         break;
		default : assert(0 && "wrong color format" ); break;
	}

	switch(dfmt)
	{
		case DEPTH16 : mDepth_internalFormat = GL_DEPTH_COMPONENT16; mDepth_format = GL_DEPTH_COMPONENT; mDepth_type = GL_FLOAT; break;
		case DEPTH24 : mDepth_internalFormat = GL_DEPTH_COMPONENT24; mDepth_format = GL_DEPTH_COMPONENT; mDepth_type = GL_FLOAT; break;
		case DEPTH32 : mDepth_internalFormat = GL_DEPTH_COMPONENT32; mDepth_format = GL_DEPTH_COMPONENT; mDepth_type = GL_FLOAT; break;
		case DEPTH24STENCIL8 : mDepth_internalFormat = GL_DEPTH24_STENCIL8; mDepth_format = GL_DEPTH_STENCIL; mDepth_type = GL_UNSIGNED_INT_24_8; break;
		default : assert(0 && "wrong depth format" ); break;
	}


	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);


	glBindTexture(GL_TEXTURE_2D, mColorTex);
	glTexImage2D(GL_TEXTURE_2D, 0, mColor_internalFormat, maxW, maxH, 0, mColor_format, mColor_type, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);


	glBindTexture(GL_TEXTURE_2D, mDepthTex);
	glTexImage2D(GL_TEXTURE_2D, 0, mDepth_internalFormat, maxW, maxH, 0, mDepth_format, mDepth_type, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);


	glBindFramebuffer(GL_FRAMEBUFFER, mGLFBOId);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, mColorTex, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  GL_TEXTURE_2D, mDepthTex, 0);
	if( dfmt == DEPTH24STENCIL8 ) {
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_STENCIL_ATTACHMENT,  GL_TEXTURE_2D, mDepthTex, 0);
	}
	checkGLFrameBufferStatus();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// NECESSARY..... else gl will align it on 4 and expect bigger buffers when readpixels is called
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	mMaxWidth  = maxW;
	mMaxHeight = maxH;

	checkGLErrors("OFFSCREENRENDERER : end init");
}




SilGLRenderJob::SilGLRenderJob(int w, int h, const double* GLPMat, const double* GLMMat, const double* GLTMat, GLfloat* vBuffer, GLuint* iBuffer, int iBuffersize)
{
	mW = w;
	mH = h;
	mVBuffer = vBuffer;
	mIBuffer = iBuffer;
	mIBufferSize = iBuffersize;
	std::copy(GLPMat, GLPMat + 16, mGLPMat);
	std::copy(GLMMat, GLMMat + 16, mGLMMat);
	std::copy(GLTMat, GLTMat + 16, mGLTMat);
}

void SilGLRenderJob::run() const
{
	glViewport(0, 0, mW, mH);

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
#ifndef DEBUG_OPENGL
	glLoadMatrixd(mGLPMat);
#endif // !DEBUG_OPENGL		

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

#ifndef DEBUG_OPENGL
	glLoadMatrixd(mGLMMat);
#else
	glTranslatef(-0.5, -0.5, -2.0);
#endif // !DEBUG_OPENGL	

	glEnable(GL_DEPTH_TEST);
	glEnableClientState(GL_VERTEX_ARRAY);
#ifdef DEBUG_OPENGL

	GLfloat vbuffer[] = {
		0, 0, 0,
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
		1, 1, 0,
		1, 0, 1,
		0, 1, 1,
		1, 1, 1
	};
	GLuint ibuffer[] = {
		0, 1, 4, 2,
		0, 2, 6, 3,
		0, 3, 5, 1,
		1, 4, 7, 5,
		5, 7, 6, 3,
		2, 6, 7, 4
	};
	GLfloat cbuffer[] = {
		0, 0, 0,
		255, 255, 255,
		0, 0, 0,
		255, 255, 255,
		255, 255, 255,
		255, 255, 255,
		0, 0, 0,
		255, 255, 255
	};
	glColorPointer(3, GL_FLOAT, 0, cbuffer);
	glVertexPointer(3, GL_FLOAT, 0, vbuffer);
	glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, ibuffer);
#else

	glVertexPointer(3, GL_FLOAT, 0, mVBuffer);
	glDrawElements(GL_TRIANGLES, mIBufferSize, GL_UNSIGNED_INT, mIBuffer);
#endif // DEBUG_OPENGL

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisable(GL_VERTEX_ARRAY);
	glDisable(GL_DEPTH_TEST);
}


ColorGLRenderJob::ColorGLRenderJob(int w, int h,
	const double* GLPMat, const double* GLMMat, const double* GLTMat,
	GLfloat* vBuffer, GLfloat* cBuffer,
	GLuint* iBuffer, int iBuffersize)
{
	mW = w;
	mH = h;
	mVBuffer = vBuffer;
	mIBuffer = iBuffer;
	mCBuffer = cBuffer;
	mIBufferSize = iBuffersize;
	std::copy(GLPMat, GLPMat + 16, mGLPMat);
	std::copy(GLMMat, GLMMat + 16, mGLMMat);
	std::copy(GLTMat, GLTMat + 16, mGLTMat);
}

void ColorGLRenderJob::run() const
{
	glViewport(0, 0, mW, mH);

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
#ifndef DEBUG_OPENGL
	glLoadMatrixd(mGLPMat);
#endif // !DEBUG_OPENGL		

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

#ifndef DEBUG_OPENGL
	glLoadMatrixd(mGLMMat);
#else
	glTranslatef(-0.5, -0.5, -2.0);
#endif // !DEBUG_OPENGL	

	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_VERTEX_ARRAY);		
	//glVertexPointer(3,GL_FLOAT, 3*sizeof(GLfloat),mVBuffer);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glColor3f(1, 1, 1);
#ifdef DEBUG_OPENGL

	GLfloat vbuffer[] = {
		0, 0, 0,
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
		1, 1, 0,
		1, 0, 1,
		0, 1, 1,
		1, 1, 1
	};
	GLuint ibuffer[] = {
		0, 1, 4, 2,
		0, 2, 6, 3,
		0, 3, 5, 1,
		1, 4, 7, 5,
		5, 7, 6, 3,
		2, 6, 7, 4
	};
	GLfloat cbuffer[] = {
		0, 0, 0,
		255, 255, 255,
		0, 0, 0,
		255, 255, 255,
		255, 255, 255,
		255, 255, 255,
		0, 0, 0,
		255, 255, 255
	};
	glColorPointer(3, GL_FLOAT, 0, cbuffer);
	glVertexPointer(3, GL_FLOAT, 0, vbuffer);
	glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, ibuffer);
#else

	glColorPointer(3, GL_FLOAT, 0, mCBuffer);
	glVertexPointer(3, GL_FLOAT, 0, mVBuffer);
	glDrawElements(GL_TRIANGLES, mIBufferSize, GL_UNSIGNED_INT, mIBuffer);
#endif // DEBUG_OPENGL


	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisable(GL_VERTEX_ARRAY);

	glDisable(GL_DEPTH_TEST);
}



}
