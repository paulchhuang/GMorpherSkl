/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <MVLib.h>
#include <time.h>
#include <iostream>

namespace MVLib
{

	/**
	 * Auxilliary function...
	 * takes a binary image as input, and computes the places where we have edges
	 */
	void findEdges(const int w, const int h, unsigned char* const src, unsigned char* dst)
	{
		int ul = -w - 1; int uc = -w; int ur = -w + 1;
		int cl = - 1;                 int cr = -w + 1;
		int ll =  w - 1; int lc = w;  int lr =  w + 1;

		unsigned char *src_ptr = src + w + 1;
		unsigned char *dst_ptr = dst + w + 1;

		std::fill(dst, dst + w*h, 0x00);
		int xmin = 1;
		int ymin = 1;
		int xmax = w-1;
		int ymax = h-1;

		
		for(int y = ymin; y < ymax; ++y){
			for(int x = xmin; x < xmax; ++x){
				if(*src_ptr){ // if we are inside the silhouette
					*dst_ptr =  ~ (src_ptr[ul] & src_ptr[uc] & src_ptr[ur] &
								   src_ptr[cl] &               src_ptr[cr] &
								   src_ptr[ll] & src_ptr[lc] & src_ptr[lr] );
				}
				src_ptr++;
				dst_ptr++;
			}
			src_ptr += 2;
			dst_ptr += 2;
		}
		
	}

	void findEdges(const int w, const int h, unsigned char* const src, bool* dst)
	{
		int ul = -w - 1; int uc = -w; int ur = -w + 1;
		int cl = - 1;                 int cr = -w + 1;
		int ll =  w - 1; int lc = w;  int lr =  w + 1;

		std::fill(dst, dst + w*h, false);

		unsigned char *src_ptr = src + w + 1;
		bool *dst_ptr = dst + w + 1;

		
		int xmin = 1;
		int ymin = 1;
		int xmax = w-1;
		int ymax = h-1;

		
		for(int y = ymin; y < ymax; ++y){
			for(int x = xmin; x < xmax; ++x){
				if(*src_ptr){ // true, as long as *src_ptr is not 0x00
					*dst_ptr =  ! (((int)src_ptr[ul]>0) & ((int)src_ptr[uc]>0) & ((int)src_ptr[ur]>0) &
								   ((int)src_ptr[cl]>0) &						 ((int)src_ptr[cr]>0) &
								   ((int)src_ptr[ll]>0) & ((int)src_ptr[lc]>0) & ((int)src_ptr[lr]>0) );
				}
				src_ptr++;
				dst_ptr++;
			}
			src_ptr += 2;
			dst_ptr += 2;
		}
		
	}

}