/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <PatchedCloud.h>

namespace SklRgedPatchedCloud 
{

using CC3D::make_float3;

// ######################################################
// ######################################################
// INTERPOLATING
// ######################################################
// ######################################################
void build_w_smooth( int                     numPatches,
                     const double            distVariance,
                     const int*              sadj_bounds,
                     const int*              pv_bounds,
                     const int*              smooth_bounds,
                     const float3*           dX0_smooth,
                     float*                  w_smooth )
{
	int numVertices_smooth = smooth_bounds[numPatches];
	int numVertices        = pv_bounds[numPatches];


	// -------------------------------------
	// 0 - allocate sum buffer
	float* sumWeights = new float[numVertices];
	std::fill(sumWeights, sumWeights+numVertices, 0.);

	// -------------------------------------
	// 1 - write unnormalised weights
	for(int pi=0;pi<numPatches;++pi)
	{
		float* w_itr        = w_smooth + smooth_bounds[pi];
		int numNeighbours_i = sadj_bounds[pi+1] - sadj_bounds[pi];
		int numVertices_i   = pv_bounds[pi+1] - pv_bounds[pi];
		const float3* const dX0_end = dX0_smooth + (numNeighbours_i+1)*numVertices_i;

		while( dX0_smooth != dX0_end ) {
			float norm = dot( *dX0_smooth, *dX0_smooth);
			*w_itr = exp( - norm / distVariance );
			w_itr       ++;
			dX0_smooth ++;
		}
	}

	// -------------------------------------
	// 2 - accum weights
	for(int pi = 0; pi < numPatches; ++pi)
	{
		float* w_itr        = w_smooth + smooth_bounds[pi];
		int numNeighbours_i = sadj_bounds[pi+1] - sadj_bounds[pi];
		float* w_sum_begin  = sumWeights + pv_bounds[pi];
		float* w_sum_end    = sumWeights + pv_bounds[pi+1];
		// sum our own
		{
			float* w_sum_itr = w_sum_begin;
			while( w_sum_itr != w_sum_end )  { *w_sum_itr += *w_itr; w_sum_itr++; w_itr++; }
		}
		// sum others
		for(int nj=0;nj<numNeighbours_i;++nj)
		{
			float* w_sum_itr = w_sum_begin;
			while( w_sum_itr != w_sum_end )  { *w_sum_itr += *w_itr; w_sum_itr++; w_itr++; }
		}
	}

	// -------------------------------------
	// 3 - normalize weights
	for(int pi = 0; pi < numPatches; ++pi)
	{
		float* w_itr        = w_smooth + smooth_bounds[pi];
		int numNeighbours_i = sadj_bounds[pi+1] - sadj_bounds[pi];
		float* w_sum_begin  = sumWeights + pv_bounds[pi];
		float* w_sum_end    = sumWeights + pv_bounds[pi+1];
		// sum our own
		{
			float* w_sum_itr = w_sum_begin;
			while( w_sum_itr != w_sum_end ){  *w_itr /= *w_sum_itr; w_itr++; w_sum_itr++; }
		}
		// sum others
		for(int nj=0;nj<numNeighbours_i;++nj)
		{
			float* w_sum_itr = w_sum_begin;
			while( w_sum_itr != w_sum_end ) {  *w_itr /= *w_sum_itr; w_itr++; w_sum_itr++; }
		}
	}


	delete[] sumWeights;
}








void interpolate_cloud_smooth( int           numPatches,
                               const float*  RT,
                               const int*    sadj,
                               const int*    sadj_bounds,
                               const int*    pv_bounds,
                               const float3* dX0_smooth,
                               const float*  w_smooth,
                               float3*       X )
{
	int numVertices = pv_bounds[numPatches];

	// reset the output vector
	float3 zero = make_float3(0,0,0);
	std::fill( X, X + numVertices, zero );

	for(int pi=0;pi<numPatches;++pi ) 
	{
		const int* nBegin = sadj + sadj_bounds[pi];
		const int* nEnd   = sadj + sadj_bounds[pi+1];
		int numVertices_i = pv_bounds[pi+1] - pv_bounds[pi];
		float3* const Xi_begin = X + pv_bounds[pi];
		// us
		{
			const float3* const dX0_end = dX0_smooth + numVertices_i;
			const float* Ri  = RT + 12*pi;
			const float3& ti = *(const float3*)(Ri+9);
			float3* X_itr = Xi_begin;
			while( dX0_smooth != dX0_end ) *X_itr++ +=  (*w_smooth++) * ( prod3(Ri, *dX0_smooth++) + ti );
//			while( dX0_smooth != dX0_end ) *X_itr++ +=  ( prod3(Ri, *dX0_smooth++) + ti );

		}
		// us
		for(const int* n_itr = nBegin; n_itr != nEnd; ++n_itr )
		{
			const float3* const dX0_end = dX0_smooth + numVertices_i;
			const float* Rj  = RT + 12*(*n_itr);
			const float3& tj = *(const float3*)(Rj+9);
			float3* X_itr = Xi_begin;
			while( dX0_smooth != dX0_end ) *X_itr++ += (*w_smooth++) * ( prod3(Rj, *dX0_smooth++) + tj );
//			dX0_smooth = dX0_end;
		}
	}
}



} // end namespace PatchedCloud


//void interpolate_Cloudsmooth( int           numPatches,
//                              const int*    sadj_bounds, // know how many neighbours
//                              const int*    pv_bounds,         // know how many vertices
//                              const float3* XN_smooth,
//                              const float*  w_smooth,
//                              float3*       X)
//{
//	// -----------
//	// 1 - reset the vector
//	float3 zero = make_float3(0,0,0);
//	std::fill( X, X + pv_bounds[numPatches], zero );

//	// -----------
//	// 2 - accum
//	for( int pi=0;pi<numPatches;++pi)
//	{
//		int numNeighbours_i = sadj_bounds[pi+1] - sadj_bounds[pi];
//		int numVertices_i   = pv_bounds[pi+1] - pv_bounds[pi];

//		float3* const Xi_end = X + numVertices_i;
//		// for us
//		{
//			float3* Xi_itr = X;
//			while (Xi_itr != Xi_end ){
//				*Xi_itr += *w_smooth * *XN_smooth;
//				Xi_itr    ++;
//				XN_smooth +=2;
//				w_smooth  ++;
//			}
//		}

//		// for the neighbours
//		for(int nj=0;nj<numNeighbours_i;++nj)
//		{
//			float3* Xi_itr = X;
//			while (Xi_itr != Xi_end ){
//				*Xi_itr += *w_smooth * *XN_smooth;
//				Xi_itr    ++;
//				XN_smooth +=2;
//				w_smooth  ++;
//			}
//		}

//		X += numVertices_i;
//	}
//}

