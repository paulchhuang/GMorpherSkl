/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <MVLib.h>
#include <cstdio>

namespace MVLib
{

	std::string buildFilename(const char* baseName, int id)
	{
		char cfilename[1024]; //buffer for the filename
		sprintf(cfilename, baseName, id);
		return std::string(cfilename);
	}

} // end namespace MVLib


