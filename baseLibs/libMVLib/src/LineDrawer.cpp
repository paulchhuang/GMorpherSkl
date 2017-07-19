/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef LINERASTERIZER_H_DEFINED
#define LINERASTERIZER_H_DEFINED

#include <cassert>
#include <cmath>
#include <cstdlib>

namespace MVLib
{

template <class T>
class LineRasterizerDXsupDY
{
	protected :
	int dx,dy;
	T*  _image;
	int _linewidth;
	T*   _Curr;
	int  _x_inc;
	int  _y_inc;
	int _a;

	public :
	inline LineRasterizerDXsupDY(T* image, int x1, int y1, int x2, int y2, int lineWidth, int deltaX)
	{
		_image = image;
		_linewidth = lineWidth;
		_Curr = image + deltaX*x1 + lineWidth*y1;

		dx = x2-x1;
		if(dx<0) dx = -dx;
		dx++;
		dy = y2-y1;
		if(dy<0) dy = -dy;
		dy++;

		assert( dx >= dy);

		_a = dx;
		_a -= dy/2; /* better symmetry */
		_y_inc = (y2 >= y1) ? lineWidth : -lineWidth;
		_x_inc = (x2 >= x1) ? deltaX : -deltaX;
	}

	inline int offset() { return _Curr - _image; }
	inline T& operator * () { return *_Curr; }
	inline T* operator () () { return _Curr; }

	inline void operator -- ()
	{
		_Curr -= _x_inc;
		_a -=dy;
		if(_a < 0)
		{
			_a += dx;
			_Curr -= _y_inc;
		}
	}

	inline void operator ++ ()
	{
		_Curr += _x_inc;
		_a -=dy;
		if(_a < 0)
		{
			_a += dx;
			_Curr += _y_inc;
		}
	}
};



template <class T>
class LineRasterizerDYsupDX
{
	protected :
	int dx,dy;
	T*  _image;
	int _linewidth;
	T*   _Curr;
	int  _x_inc;
	int  _y_inc;
	int _a;

	public :
	inline LineRasterizerDYsupDX(T* image, int x1, int y1, int x2, int y2, int lineWidth, int deltaX)
	{
		_image = image;
		_linewidth = lineWidth;
		_Curr = image + deltaX*x1 + lineWidth*y1;

		dx = x2-x1;
		if(dx<0) dx = -dx;
		dx++;
		dy = y2-y1;
		if(dy<0) dy = -dy;
		dy++;

		assert( dy >= dx);

		_a = dy;
		_a -= dx/2; /* better symmetry */
		_y_inc = (y2 >= y1) ? lineWidth : -lineWidth;
		_x_inc = (x2 >= x1) ? deltaX : -deltaX;
	}

	inline int offset() { return _Curr - _image; }
	inline T& operator * () { return *_Curr; }
	inline T* operator () (){ return _Curr; }

	inline void operator -- ()
	{
		_Curr -= _y_inc;
		_a -=dx;
		if(_a < 0)
		{
			_a += dy;
			_Curr -= _x_inc;
		}
	}

	inline void operator ++ ()
	{
		_Curr += _y_inc;
		_a -= dx;
		if(_a < 0)
		{
			_a += dy;
			_Curr += _x_inc;
		}
	}
};



void drawLineGRAY( unsigned char* image, int W, int x1, int y1, int x2, int y2, unsigned char val )
{
	if( abs(x1-x2) >= abs(y2-y1) )
	{
		LineRasterizerDXsupDY<unsigned char> r(image , x1, y1, x2, y2, W, 1);
		int n = abs(x2 - x1);
		for(int i=0;i<n;++i){
			*r = 1;
			++r;}
	}
	else
	{
		LineRasterizerDYsupDX<unsigned char> r(image , x1, y1, x2, y2, W, 1);
		int n = abs(y2-y1);
		for(int i=0;i<n;++i){
			*r = 1;
			++r;}
	}
}


void drawLineRGB( unsigned char* image, int W, int x1, int y1, int x2, int y2, unsigned char r, unsigned char g, unsigned char b )
{
	if( abs(x1-x2) >= abs(y2-y1) )
	{
		LineRasterizerDXsupDY<unsigned char> ritr(image , x1, y1, x2, y2, W*3, 3);
		int n = abs(x2 - x1);
		for(int i=0;i<n;++i){
			ritr()[0] = r;
			ritr()[1] = g;
			ritr()[2] = b;
			++ritr;}
	}
	else
	{
		LineRasterizerDYsupDX<unsigned char> ritr(image , x1, y1, x2, y2, W*3, 3);
		int n = abs(y2-y1);
		for(int i=0;i<n;++i){
			ritr()[0] = r;
			ritr()[1] = g;
			ritr()[2] = b;
			++ritr;}
	}
}


} // end namespace MVLib


#endif
