/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <PatchedCloud.h>


//TO DO: make all this transform run in parallel


namespace SklRgedPatchedCloud 
{


//void transform_Cloud( int           numPatches,
//                      const float*  RT,
//                      const int*    pv_bounds,
//                      const float3* dXN0,
//                      float3*       XN )
//{
//	float3* XN_itr = (float3*)XN;  
//	for(int pi=0;pi<numPatches;++pi) {
//		const float3*       dXN0_itr =  dXN0 + 2*pv_bounds[pi];
//		const float3* const dXN0_end =  dXN0 + 2*pv_bounds[pi+1];
//
//		const float* R = RT + 12*pi;
//		const float3& t = *(const float3*)(R + 9);
//		while ( dXN0_itr != dXN0_end ) {
//			// coord
//			XN_itr[0] = prod3( R, dXN0_itr[0] ) + t; // coord
//			XN_itr[1] = prod3( R, dXN0_itr[1] ); // normal
//			XN_itr   +=2;
//			dXN0_itr +=2;
//		}
//	}
//}

// parallel version
void transform_Cloud( int           numPatches,
                      const float*  RT,
                      const int*    pv_bounds,
                      const float3* dXN0,
                      float3*       XN )
{
	#pragma omp parallel for 
	for(int pi=0;pi<numPatches;++pi) {
		
		const float3*       dXN0_itr =  dXN0 + 2*pv_bounds[pi];
		const float3* const dXN0_end =  dXN0 + 2*pv_bounds[pi+1];
		float3*				XN_itr	 =  XN + 2*pv_bounds[pi];  

		const float* R = RT + 12*pi;
		const float3& t = *(const float3*)(R + 9);
		while ( dXN0_itr != dXN0_end ) {
			// coord
			XN_itr[0] = prod3( R, dXN0_itr[0] ) + t; // coord
			XN_itr[1] = prod3( R, dXN0_itr[1] ); // normal
			XN_itr   +=2;
			dXN0_itr +=2;
		}

		assert(XN_itr==(XN+2*pv_bounds[pi+1]));
	}
}


//void rotate_Cloud( int           numPatches,
//                   const float*  RT,
//                   const int*    pv_bounds,
//                   const float3* dXN0,
//                   float3*       dXN )
//{
//	float3* dXN_itr = (float3*)dXN;	
//	for(int pi=0;pi<numPatches;++pi) {
//		const float3*       dXN0_itr =  dXN0 + 2*pv_bounds[pi];
//		const float3* const dXN0_end =  dXN0 + 2*pv_bounds[pi+1];
//
//		const float* R = RT + 12*pi;
//		while ( dXN0_itr != dXN0_end ) {
//			dXN_itr[0] = prod3( R, dXN0_itr[0] ); // coord
//			dXN_itr[1] = prod3( R, dXN0_itr[1] ); // normal
//			dXN_itr   +=2;
//			dXN0_itr +=2;
//		}
//	}
//}

// parallel version
void rotate_Cloud( int           numPatches,
                   const float*  RT,
                   const int*    pv_bounds,
                   const float3* dXN0,
                   float3*       dXN )
{
	#pragma omp parallel for 
	for(int pi=0;pi<numPatches;++pi) {
		const float3*       dXN0_itr =  dXN0 + 2*pv_bounds[pi];
		const float3* const dXN0_end =  dXN0 + 2*pv_bounds[pi+1];
		float3*				dXN_itr  =  dXN + 2*pv_bounds[pi];	
		const float* R = RT + 12*pi;
		while ( dXN0_itr != dXN0_end ) {
			dXN_itr[0] = prod3( R, dXN0_itr[0] ); // coord
			dXN_itr[1] = prod3( R, dXN0_itr[1] ); // normal
			dXN_itr   +=2;
			dXN0_itr +=2;
		}
		assert(dXN_itr==(dXN + 2*pv_bounds[pi+1]));
	}
}

// ######################################################
// SMOOTH
void transform_Cloudsmooth( int           numPatches,
                           const float*  RT,
                           const int*    sadj,
                           const int*    sadj_bounds,
                           const int*    pv_bounds,
                           const float3* dXN0_smooth,
                           float3*       XN_smooth )
{
	float3* XN_itr         = XN_smooth;
	const float3* dXN0_itr = dXN0_smooth;	
	for(int pi=0;pi<numPatches;++pi)
	{
		int numVertices_i = pv_bounds[pi+1] - pv_bounds[pi];
		const int*    const nBegin    = sadj + sadj_bounds[pi];
		const int*    const nEnd      = sadj + sadj_bounds[pi+1];
		// ------------------
		// FOR US
		{
			const float3* const dXN0_end  = dXN0_itr + 2*numVertices_i;
			const float* R = RT + 12*pi;
			const float3& t = *(const float3*)(R + 9);			
			while ( dXN0_itr != dXN0_end ) {
				XN_itr[0] = prod3(R, dXN0_itr[0] ) + t; // coord
				XN_itr[1] = prod3(R, dXN0_itr[1] );     // normal
				XN_itr   +=2;
				dXN0_itr +=2;				
			}
		}
		// ------------------
		// FOR NEIGHBOURS
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3* const dXN0_end  = dXN0_itr + 2*numVertices_i;
			const float* R = RT + 12*(*pj_itr);
			const float3& t = *(const float3*)(R + 9);
			while ( dXN0_itr != dXN0_end ) {
				XN_itr[0] = prod3(R, dXN0_itr[0] ) + t; // coord
				XN_itr[1] = prod3(R, dXN0_itr[1] );     // normal
				XN_itr   +=2;
				dXN0_itr +=2;
			}
		}
	}
}
// parallel version
void transform_Cloudsmooth( int           numPatches,
                           const float*  RT,
						   const int*	smooth_bound,
                           const int*    sadj,
                           const int*    sadj_bounds,
                           const int*    pv_bounds,
                           const float3* dXN0_smooth,
						   float3*       XN_smooth ){
	
	#pragma omp parallel for 
	for(int pi=0;pi<numPatches;++pi)
	{
		int numVertices_i = pv_bounds[pi+1] - pv_bounds[pi];
		const int*    const nBegin    = sadj + sadj_bounds[pi];
		const int*    const nEnd      = sadj + sadj_bounds[pi+1];

		const float3* dXN0_itr = dXN0_smooth + 2*smooth_bound[pi];
		float3* XN_itr         = XN_smooth + 2*smooth_bound[pi];			

		// ------------------
		// FOR US
		{
			const float3* const dXN0_end  = dXN0_itr + 2*numVertices_i;
			const float* R = RT + 12*pi;
			const float3& t = *(const float3*)(R + 9);			
			while ( dXN0_itr != dXN0_end ) {
				XN_itr[0] = prod3(R, dXN0_itr[0] ) + t; // coord
				XN_itr[1] = prod3(R, dXN0_itr[1] );     // normal
				XN_itr   +=2;
				dXN0_itr +=2;				
			}
		}
		// ------------------
		// FOR NEIGHBOURS
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3* const dXN0_end  = dXN0_itr + 2*numVertices_i;
			const float* R = RT + 12*(*pj_itr);
			const float3& t = *(const float3*)(R + 9);
			while ( dXN0_itr != dXN0_end ) {
				XN_itr[0] = prod3(R, dXN0_itr[0] ) + t; // coord
				XN_itr[1] = prod3(R, dXN0_itr[1] );     // normal
				XN_itr   +=2;
				dXN0_itr +=2;
			}
		}

		assert(dXN0_itr==(dXN0_smooth + 2*smooth_bound[pi+1]));
		assert(XN_itr==(XN_smooth + 2*smooth_bound[pi+1]));
	}
}


void rotate_dX_smooth( int                     numPatches,
                       const float*            RTf,
                       const int*              sadj,
                       const int*              sadj_bounds,
                       const int*              pv_bounds,
                       const float3*           dX0_smooth,
                       float3*                 dX_smooth )
{
	for(int pi=0;pi<numPatches;++pi)
	{
		int numVertices_i = pv_bounds[pi+1] - pv_bounds[pi];
		const int* nBegin = sadj + sadj_bounds[pi];
		const int* nEnd   = sadj + sadj_bounds[pi+1];
		// ------------------
		// FOR US
		{
			const float3* const dX0_end  = dX0_smooth + numVertices_i;
			const float* R = RTf + 12*pi;
			while ( dX0_smooth != dX0_end ) *dX_smooth++ = prod3( R, *dX0_smooth++);
		}
		// ------------------
		// FOR NEIGHBOURS
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3* const dX0_end  = dX0_smooth + numVertices_i;
			const float* R = RTf + 12*(*pj_itr);
			while ( dX0_smooth != dX0_end ) *dX_smooth++ = prod3( R, *dX0_smooth++);
		}
	}
}
// parallel version
void rotate_dX_smooth( int                     numPatches,					   
                       const float*            RTf,
					   const int*			   smooth_bound,
                       const int*              sadj,
                       const int*              sadj_bounds,
                       const int*              pv_bounds,
                       const float3*           dX0_smooth,
                       float3*                 dX_smooth )
{

	#pragma omp parallel for 
	for(int pi=0;pi<numPatches;++pi)
	{
		int numVertices_i = pv_bounds[pi+1] - pv_bounds[pi];
		const int* nBegin = sadj + sadj_bounds[pi];
		const int* nEnd   = sadj + sadj_bounds[pi+1];

		const float3*           dX0_smooth_pi	= dX0_smooth + smooth_bound[pi];
        float3*                 dX_smooth_pi	= dX_smooth + smooth_bound[pi];
		// ------------------
		// FOR US
		{
			const float3* const dX0_end  = dX0_smooth_pi + numVertices_i;
			const float* R = RTf + 12*pi;
			while ( dX0_smooth_pi != dX0_end ) *dX_smooth_pi++ = prod3( R, *dX0_smooth_pi++);
		}
		// ------------------
		// FOR NEIGHBOURS
		for(const int* pj_itr = nBegin; pj_itr != nEnd; ++pj_itr )
		{
			const float3* const dX0_end  = dX0_smooth_pi + numVertices_i;
			const float* R = RTf + 12*(*pj_itr);
			while ( dX0_smooth_pi != dX0_end ) *dX_smooth_pi++ = prod3( R, *dX0_smooth_pi++);
		}

		assert(dX0_smooth_pi==(dX0_smooth + smooth_bound[pi+1]));
		assert(dX_smooth_pi==(dX_smooth + smooth_bound[pi+1]));
	}
}


} // end namespace PatchedCloud
