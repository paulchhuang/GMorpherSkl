/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef KERNEL_ICP_H_DEFINED
#define KERNEL_ICP_H_DEFINED


#define CUDA_ICP_ESTEP_TPB 32


void runKernel_ICP_EStep( int         numPatches,
                          int         numTargets,
                          float       sigma,
                          float       dotThresh,
                          float       EOutlier,
                          const void* XNt_device,
                          const void* XN_device,
                          const void* XN_bounds_device,
                          void*       weights_device,
                          void*       targets_device );


#endif
