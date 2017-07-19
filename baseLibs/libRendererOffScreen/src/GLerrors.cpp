/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/

#include <windows.h>

//#include <GL/gl.h>
//#include <GL/glu.h>
#include <GL/glew.h>
//#include <GL/glext.h>
#include "GLerrors.h"
#include <iostream>

namespace OffScreenRendering
{



	// #########################################################
	// #########################################################
	// ERROR CHECKING STUFF
	// #########################################################
	// #########################################################


	void checkGLErrors(const char *label)
	{
		GLenum errCode;
		const GLubyte *errStr;
		if ((errCode = glGetError()) != GL_NO_ERROR)
		{
			errStr = gluErrorString(errCode);
			std::cout<<"OpenGL ERROR: ";
			std::cout<<errStr;
			std::cout<<"Label: ";
			std::cout<<label<<std::endl;
		}
	}





	void checkGLFrameBufferStatus()
	{
		GLenum status;
		status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		switch (status)
		{
			case GL_FRAMEBUFFER_COMPLETE:
				//std::cout<<"GL_FRAMEBUFFER_COMPLETE_EXT"<<std::endl;
				break;
			case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
				std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT"<<std::endl;
				break;
			case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
				std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT"<<std::endl;
				break;
			//case GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT:
			//	std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT"<<std::endl;
			//	break;
			//~ case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS:
				//~ std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS"<<std::endl;
				//~ break;
			//~ case GL_FRAMEBUFFER_INCOMPLETE_FORMATS:
				//~ std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_FORMATS"<<std::endl;
				//~ break;
			case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
				std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER"<<std::endl;
				break;
			case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
				std::cout<<"GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER"<<std::endl;
				break;
			case GL_FRAMEBUFFER_UNSUPPORTED:
				std::cout<<"GL_FRAMEBUFFER_UNSUPPORTED"<<std::endl;
				break;
			//case GL_FRAMEBUFFER_STATUS_ERROR_EXT:
			//	std::cout<<"GL_FRAMEBUFFER_STATUS_ERROR_EXT"<<std::endl;
			//	break;
			default:
				std::cout<<"UNKNOWN FBO ERROR"<<std::endl;
				break;
		}
	}


} // end namespace OffScreenRendering


