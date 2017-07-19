/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef LIBOFFSCREENRENDERING_IGLRENDERJOB_H_DEFINED
#define LIBOFFSCREENRENDERING_IGLRENDERJOB_H_DEFINED


namespace OffScreenRendering
{
	class IGLRenderJob
	{
		public :
		virtual ~IGLRenderJob(){};
		virtual void run() const = 0;
	};
}



#endif
