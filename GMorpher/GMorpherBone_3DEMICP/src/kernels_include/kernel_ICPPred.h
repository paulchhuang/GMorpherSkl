/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef KERNEL_ICPPRED_H_DEFINED
#define KERNEL_ICPPRED_H_DEFINED


#define CUDA_ICPPRED_ESTEP_TPB 64

// careful this function returns indices in the smooth vector...
// there needs to be some post-processing on them to map to the normal vector
void runKernel_ICPPred_EStep( int         numPatches,
                              int         numTargets,
                              float       sigma,
                              float       dotThresh,
                              float       EOutlier,
                              const void* XNt_device,
                              const void* XN_smooth_device,
                              const void* XN_smooth_bounds_device,
                              void*       weights_device,
                              void*       targets_device);

#endif
