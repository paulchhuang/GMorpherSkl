/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef LIBOFFSCREENRENDERING_RENDERER_H_DEFINED
#define LIBOFFSCREENRENDERING_RENDERER_H_DEFINED

#include <windows.h> 
#include <gl\glew.h>
#include "IGLRenderJob.h"


namespace OffScreenRendering
{

	class Renderer
	{
		public :
		enum ColorFormats{ RGB8, RGBA8, GRAY8, GRAY32};
		enum DepthFormats{ DEPTH16, DEPTH24, DEPTH32, DEPTH24STENCIL8};

		Renderer(unsigned int maxW, unsigned int maxH,
		         ColorFormats cfmt, DepthFormats dfmt);

		~Renderer();

		void resize(unsigned int maxW, unsigned int maxH,
		            ColorFormats cfmt, DepthFormats dfmt);

		void draw(IGLRenderJob* job,
		          unsigned int x, unsigned int y, unsigned int dx, unsigned int dy,
		          unsigned char* out_image,
		          float*         out_depth);

		void draw(IGLRenderJob* job,
		          unsigned int x, unsigned int y, unsigned int dx, unsigned int dy,
		          float*         out_image,
		          float*         out_depth);

		void draw(IGLRenderJob* job,
		          unsigned int x, unsigned int y, unsigned int dx, unsigned int dy,
		          unsigned char* out_image,
		          float*         out_depth,
		          unsigned char* out_stencil);


		protected :

		void allocate(unsigned int maxW, unsigned int maxH,
		              ColorFormats cfmt, DepthFormats dfmt);

		unsigned int           mMaxWidth,mMaxHeight; //holds the image size
		GLuint        mGLFBOId;
		GLuint        mColorTex;
		GLuint        mDepthTex;

		GLint         mColor_internalFormat;
		GLint         mColor_format;
		GLint         mColor_type;

		GLint         mDepth_internalFormat;
		GLint         mDepth_format;
		GLint         mDepth_type;
	};


	/**
	* Auxilliary class...
	* is a simple renderer
	*/
	class SilGLRenderJob : public IGLRenderJob
	{ // this class deals with BW silhouette rendering (no color information). For color rendering see energy_imgblob
		protected:
			double  mGLPMat[16], mGLMMat[16], mGLTMat[16];
			int mW, mH;
			GLfloat* mVBuffer;
			GLuint* mIBuffer;
			int mIBufferSize;

		public:
			SilGLRenderJob(int w, int h, const double* GLPMat, const double* GLMMat, const double* GLTMat, GLfloat* vBuffer, GLuint* iBuffer, int iBuffersize);
		
			virtual void run() const;			
	};


	/**
	* Auxilliary class...
	* is a simple renderer
	*/
	class ColorGLRenderJob : public IGLRenderJob
	{
	protected:
		double  mGLPMat[16], mGLMMat[16], mGLTMat[16];
		int mW, mH;
		GLfloat* mVBuffer;
		GLfloat* mCBuffer;
		GLuint* mIBuffer;
		int mIBufferSize;

	public:
		ColorGLRenderJob(int w, int h,
			const double* GLPMat, const double* GLMMat, const double* GLTMat,
			GLfloat* vBuffer, GLfloat* cBuffer,
			GLuint* iBuffer, int iBuffersize);

		virtual void run() const;
	};

	//~ class Renderer_1D
	//~ {
		//~ public :
		//~ Renderer_1D( unsigned int maxW,
		             //~ ColorFormats cfmt, DepthFormats dfmt);
		//~ ~Renderer_1D();

		//~ void draw( IGLRenderJob* job, int x, int dx,
		           //~ unsigned char* out_image,
		           //~ float*         out_depth);

		//~ void resize(unsigned int maxW,
		            //~ ColorFormats cfmt, DepthFormats dfmt);

		//~ protected :
		//~ int           mMaxWidth;
		//~ GLuint        mGLFBOId;
		//~ GLuint        mColorTex;
		//~ GLuint        mDepthTex;

		//~ GLint         mColor_internalFormat;
		//~ GLint         mColor_format;
		//~ GLint         mColor_type;

		//~ GLint         mDepth_internalFormat;
		//~ GLint         mDepth_format;
		//~ GLint         mDepth_type;
	//~ };

}

#endif
