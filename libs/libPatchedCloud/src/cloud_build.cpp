/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <PatchedCloud.h>


namespace SklRgedPatchedCloud 
{

void build_dXN0( int                                         numPatches,
                 const rigidT*                               RT,
                 const int*                                  pv_bounds,
                 const float3*                               XN0,
                 float3*                                     dXN0 )
{
	const float3* XN0_itr = XN0;
	for(int pi=0;pi<numPatches;++pi) 
	{
		const float3& ci = RT[pi].m_t;
		const float3* XN0_end = XN0 + 2*pv_bounds[pi+1];
		while( XN0_itr != XN0_end ) {
			dXN0[0] = XN0_itr[0] - ci;
			dXN0[1] = XN0_itr[1];
			XN0_itr+=2;;
			dXN0+=2;
		}
	}
}


void build_dXN0Smooth( int           numPatches,
                       const float*  RT,
                       const int*    sadj,
                       const int*    sadj_bounds,
                       const int*    pv_bounds,
                       const float3* XN0,
                       float3*       dXN0_smooth )
{
	for(int pi=0;pi<numPatches;++pi)
	{
		const int*    const nBegin    = sadj + sadj_bounds[pi];
		const int*    const nEnd      = sadj + sadj_bounds[pi+1];
		const float3* const XN0_begin  = XN0 + 2*pv_bounds[pi];
		const float3* const XN0_end    = XN0 + 2*pv_bounds[pi+1];
		const float3& ci = *( (const float3*)(RT + 12*pi+9) );
		// ------------------
		// FOR US
		{
			const float3* XN0_itr = XN0_begin;
			while ( XN0_itr != XN0_end ) {
				dXN0_smooth[0] = XN0_itr[0]  - ci; //vertex
				dXN0_smooth[1] = XN0_itr[1]; //normal
				dXN0_smooth   +=2;
				XN0_itr       +=2;
			}
		}
		// ------------------
		// FOR NEIGHBOURS
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3& cj = *( (const float3*)(RT + 12*(*pj_itr)+9) );
			const float3* XN0_itr = XN0_begin;
			while ( XN0_itr != XN0_end ) {
				dXN0_smooth[0] = XN0_itr[0]  - cj; //vertex
				dXN0_smooth[1]   = XN0_itr[1]; //normal
				dXN0_smooth   +=2;
				XN0_itr       +=2;
			}
		}
	}
}


void build_dXN0Smooth_from_dXN0( int           numPatches,
                                 const float*  RT,
                                 const int*    sadj,
                                 const int*    sadj_bounds,
                                 const int*    pv_bounds,
                                 const float3* dXN0,
                                 float3*       dXN0_smooth )
{
	float3* dXN0_smooth_itr = dXN0_smooth;

	for(int pi=0;pi<numPatches;++pi)
	{
		const int*    const nBegin    = sadj + sadj_bounds[pi];
		const int*    const nEnd      = sadj + sadj_bounds[pi+1];
		const float3* const dXN0_begin  = dXN0 + 2*pv_bounds[pi];
		const float3* const dXN0_end    = dXN0 + 2*pv_bounds[pi+1];
		const float3& ci = *( (const float3*)(RT + 12*pi+9) );
		// ------------------
		// FOR US
		{
			const float3* dXN0_itr = dXN0_begin;
			while ( dXN0_itr != dXN0_end ) {
				dXN0_smooth_itr[0] = dXN0_itr[0]; //vertex
				dXN0_smooth_itr[1] = dXN0_itr[1]; //normal
				dXN0_smooth_itr   +=2;
				dXN0_itr          +=2;
			}
		}
		// ------------------
		// FOR NEIGHBOURS
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3& cj = *( (const float3*)(RT + 12*(*pj_itr)+9) );
			float3 ci_cj = ci - cj;
			const float3* dXN0_itr = dXN0_begin;
			while ( dXN0_itr != dXN0_end ) {
				dXN0_smooth_itr[0] = dXN0_itr[0] + ci_cj; //vertex
				dXN0_smooth_itr[1] = dXN0_itr[1]; //normal
				dXN0_smooth_itr   +=2;
				dXN0_itr          +=2;
			}
		}
	}
}


void build_dX0_smooth( int                        numPatches,
                       const int*                 sadj,
                       const int*                 sadj_bounds,
                       const int*                 pv_bounds,
                       const rigidT*              RT0,
                       const float3*              X0,
                       float3*                    dX0_smooth )
{
	for( int pi=0;pi<numPatches;++pi)
	{
		const int* nBegin = sadj + sadj_bounds[pi];
		const int* nEnd   = sadj + sadj_bounds[pi+1];
		const float3* const Xi_beg = X0 + pv_bounds[pi];
		const float3* const Xi_end = X0 + pv_bounds[pi+1];

		// copy ourselves
		{
			const float3 ci = RT0[pi].m_t;
			const float3* Xi_itr = Xi_beg;
			while( Xi_itr != Xi_end )  *dX0_smooth++ = *Xi_itr++ - ci;
		}
		// copy our neigbours
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3 cj = RT0[*pj_itr].m_t;
			const float3* Xi_itr = Xi_beg;
			while(  Xi_itr != Xi_end )  *dX0_smooth++ = *Xi_itr++ - cj;
		}
	}
}


} // end namespace PatchedCloud
