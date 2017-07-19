/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <MVLib.h>


namespace MVLib
{
	void RGB_2_GRAY( const unsigned char* RGB_buffer, unsigned char* GRAY_buffer, const int w, const int h)
	{
		const int numPix = w*h;
		const int rs = 19507;
		const int gs = 38365;
		const int bs = 7153;

		const unsigned char* GRAY_end = GRAY_buffer + w*h;
		while( GRAY_buffer != GRAY_end ) {
			*GRAY_buffer = (rs*RGB_buffer[0]+gs*RGB_buffer[1]+bs*RGB_buffer[2])/65025;
			RGB_buffer+=3;
			GRAY_buffer++;
		}
	}



	void GRAY_2_RGB( const unsigned char* GRAY_buffer, unsigned char* RGB_buffer, const int w, const int h)
	{
		const unsigned char* const GRAY_end = GRAY_buffer + w*h;
		while( GRAY_buffer != GRAY_end ) {
			RGB_buffer[0] = RGB_buffer[1] = RGB_buffer[2] = *GRAY_buffer;
			RGB_buffer+=3;
			GRAY_buffer++;
		}
	}

}