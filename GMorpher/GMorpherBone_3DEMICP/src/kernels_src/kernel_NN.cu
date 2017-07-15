#include "../kernels_include/kernel_NN.h"
#include <stdio.h>
#include <float.h>
#include <cutil_math.h>



texture<int, 1, cudaReadModeElementType>    XNsmoothboundsTex;

__global__ void CUDA_NN(const int     numTargets,
						const int     numPatches,
						const float   sigma,
						const float   EOutlier,
						const float   dotThresh,
						const float3* XNt,
						const float3* XN_smooth,
						float*  weights,
						int*    targets )
{
	uint vit = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if( vit >= numTargets ) return;

	float* w_line = weights + (numPatches+1)*vit;
	int*   t_line = targets + (numPatches+1)*vit;

	float sumWeights = 0.0f;
	float  half_sigmainv2 = 0.5f/(sigma*sigma);

	const float3 Xt = XNt[2*vit];
	const float3 Nt = XNt[2*vit+1];

	const float3* xn_ptr = XN_smooth;

	w_line[numPatches] = 1;	
	for( int pi=0; pi<numPatches; ++pi)
	{
		float minVal = FLT_MAX;
		const float3* minIdx   = xn_ptr;
		const float3* const beg_ptr = xn_ptr;
		const float3* const end_ptr = XN_smooth + 2*tex1Dfetch(XNsmoothboundsTex, pi+1 );

		while( xn_ptr != end_ptr ) {
			float3 delta = xn_ptr[0] - Xt;
			float  dist2 = dot(delta,delta);
			float  dotN  = dot( xn_ptr[1], Nt );
			if( (dotN > dotThresh) && (dist2 < minVal ) ){
				minVal = dist2;
				minIdx = xn_ptr;
			}
			xn_ptr +=2;
		}

		if( minVal == FLT_MAX ) w_line[pi] = 0;
		else {
			float wi = __expf(-minVal*half_sigmainv2);
			sumWeights+= wi;
			w_line[pi] = wi;
			w_line[numPatches] *= (1 - wi);
		}
		t_line[pi] = minIdx - beg_ptr; // since we have the indices multiplied by 2
	}

	//// -----------------------------
	//// 2 - normalize
	w_line[numPatches]    = EOutlier*sigma*sigma;
	t_line[numPatches]    = -1;
	sumWeights += w_line[numPatches];

	for(int pi=0;pi<numPatches+1;++pi) w_line[pi] /= sumWeights;
}

void runKernel_ICPPred_EStep( int         numPatches,
                              int         numTargets,
                              float       sigma,
                              float       dotThresh,
                              float       EOutlier,
                              const void* XNt_device,
                              const void* XN_smooth_device,
                              const void* XN_smooth_bounds_device,
                              void*       weights_device,
                              void*       targets_device)
{
	// 1 - find the correct number of blocks
	int numBlocks   = (numTargets / CUDA_NN_TPB);
	if( numTargets % CUDA_NN_TPB ) numBlocks++;


	//printf("launching kernel with %d blocks \n", numBlocks );
	cudaEvent_t startEvt, stopEvt;
	cudaEventCreate(&startEvt);
	cudaEventCreate(&stopEvt);
	cudaEventRecord(startEvt,0);


	// ###################################
	cudaChannelFormatDesc channeldesc = cudaCreateChannelDesc<int>();
	cudaBindTexture(0, XNsmoothboundsTex, XN_smooth_bounds_device, channeldesc,  (numPatches+1)*sizeof(int));

	CUDA_PredEStep<<<numBlocks, CUDA_NN_TPB>>>( numTargets,
	                                                       numPatches,
	                                                       sigma,
	                                                       EOutlier,
	                                                       dotThresh,
	                                                       (const float3*)XNt_device,
	                                                       (const float3*)XN_smooth_device,
	                                                       (float*)       weights_device,
	                                                       (int*)         targets_device);

	
	cudaUnbindTexture( XNsmoothboundsTex);
	// ###################################

	cudaEventRecord(stopEvt,0);
	cudaThreadSynchronize();
	float ms = 0.0;
	cudaEventElapsedTime (&ms, startEvt, stopEvt);
	cudaEventDestroy(startEvt);
	cudaEventDestroy(stopEvt);

	//printf("CUDA -- spent %f ms in the kernel %s \n", ms , __FUNCTION__);

	{
		cudaError_t status = cudaGetLastError();
		if(status != cudaSuccess) printf("(E) CUDA error: %s:%d\n %s\n", __FUNCTION__, __LINE__, cudaGetErrorString(status));
	}

}
