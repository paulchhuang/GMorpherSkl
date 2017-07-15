#include "../kernels_include/kernel_ICP.h"
#include <stdio.h>
#include <float.h>
//#include <cutil_math.h>



texture<int, 1, cudaReadModeElementType>    XNboundsTex;

__global__ void CUDA_EStep(const int     numTargets,
                            const int     numPatches,
                            const float   sigma,
                            const float   EOutlier,
                            const float   dotThresh,
                            const float3* XNt,
                            const float3* XN,
                            float*  weights,
                            int*    targets )
{
	unsigned int vit = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if( vit >= numTargets ) return;

	float* w_line = weights + (numPatches+1)*vit;
	int*   t_line = targets + (numPatches+1)*vit;

	float sumWeights = 0.0f;
	float  half_sigmainv2 = 0.5f/(sigma*sigma);

	const float3 Xt = XNt[2*vit];
	const float3 Nt = XNt[2*vit+1];

	const float3* xn_ptr = XN;
	for( int pi=0; pi<numPatches; ++pi)
	{
		float minVal = FLT_MAX;
		int minIdx   = 0;
		const float3* end_ptr   = XN + 2*tex1Dfetch(XNboundsTex, pi+1 );

		while( xn_ptr != end_ptr ) {
			float3 delta = make_float3(xn_ptr[0].x - Xt.x, xn_ptr[0].y - Xt.y, xn_ptr[0].z - Xt.z);
			float  dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
			float  dotN = xn_ptr[1].x*Nt.x + xn_ptr[1].y*Nt.y + xn_ptr[1].z*Nt.z;
			if( (dotN > dotThresh) && (dist2 < minVal ) ){
				minVal = dist2;
				minIdx = xn_ptr - XN;
			}
			xn_ptr +=2;
		}

		if( minVal == FLT_MAX ) w_line[pi] = 0;
		else {
			float wi = expf(-minVal*half_sigmainv2);
			sumWeights+= wi;
			w_line[pi] = wi;
		}
		t_line[pi] = minIdx; // since we have the indices multiplied by 2
	}

	// -----------------------------
	// 2 - normalize
	w_line[numPatches]    = EOutlier*sigma*sigma;
	t_line[numPatches]    = -1;
	sumWeights += EOutlier*sigma*sigma;

	for(int pi=0;pi<numPatches+1;++pi) w_line[pi] /= sumWeights;
}




void runKernel_ICP_EStep( int         numPatches,
                          int         numTargets,
                          float       sigma,
                          float       dotThresh,
                          float       EOutlier,
                          const void* XNt_device,
                          const void* XN_device,
                          const void* XN_bounds_device,
                          void*       weights_device,
                          void*       targets_device )
{
	// 1 - find the correct number of blocks
	int numBlocks   = (numTargets / CUDA_ICP_ESTEP_TPB);
	if( numTargets % CUDA_ICP_ESTEP_TPB ) numBlocks++;


	cudaEvent_t startEvt, stopEvt;
	cudaEventCreate(&startEvt);
	cudaEventCreate(&stopEvt);
	cudaEventRecord(startEvt,0);


	// ###################################
	cudaChannelFormatDesc channeldesc = cudaCreateChannelDesc<int>();
	cudaBindTexture(0, XNboundsTex, XN_bounds_device, channeldesc,  (numPatches+1)*sizeof(int));

	CUDA_EStep<<<numBlocks, CUDA_ICP_ESTEP_TPB>>>( numTargets,
	                                               numPatches,
	                                               sigma,
	                                               EOutlier,
	                                               dotThresh,
	                                               (const float3*)XNt_device,
	                                               (const float3*)XN_device,
	                                               (float*) weights_device,
	                                               (int*) targets_device );

	cudaUnbindTexture( XNboundsTex);
	// ###################################

	cudaEventRecord(stopEvt,0);
	cudaThreadSynchronize();
	float ms = 0.0;
	cudaEventElapsedTime (&ms, startEvt, stopEvt);
	cudaEventDestroy(startEvt);
	cudaEventDestroy(stopEvt);

	printf("CUDA -- spent %f ms in the kernel %s \n", ms , __FUNCTION__);

	cudaError_t status = cudaGetLastError();
	if(status != cudaSuccess){
		printf("(E) CUDA error: %s\n %s\n", __FUNCTION__, cudaGetErrorString(status));
	}
}
