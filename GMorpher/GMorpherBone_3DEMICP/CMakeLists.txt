project (GMorpherBone_3DEMICP)
cmake_minimum_required(VERSION 2.8)


#-------------------- OpenMP part-------------------------
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)	
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# -------------------------
# THE CUDA KERNEL STUFF COMPILED WITH NVCC
set( KERNEL_SRC src/kernels_src/kernel_ICP.cu
                src/kernels_src/kernel_ICPPred.cu )
INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/include src/kernels_include)


if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
	SET( CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3   --ptxas-options=-v  -Xcompiler /MDd -D_DEBUG)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RELEASE")
	SET( CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3   --ptxas-options=-v  -Xcompiler /MD -D_RELEASE)
endif()

CUDA_COMPILE( COMPILED_KERNELS  STATIC  ${KERNEL_SRC} )

if(MSVC_VERSION STREQUAL "1800")
	set(VC_VER "12")
elseif(MSVC_VERSION STREQUAL "1900")
	set(VC_VER "14")
elseif(MSVC_VERSION STREQUAL "1910")
	set(VC_VER "15")
endif()

# -------------------------
# THE SOURCES

set( SRC_COMMON   
	 src/Solver_EM_ceres.cpp # src/Solver_EM.cpp     
	 src/IEnergyTerm_EM_ceres.cpp #src/IEnergyTerm_EM.cpp
	 src/EnergyTerm_EMPred_ceres.cpp	 #src/EnergyTerm_EM.cpp #src/EnergyTerm_EMPred.cpp
   )

INCLUDE_DIRECTORIES( include
					 # ${libKineben_SOURCE_DIR}/include
					 # ${libKinebenPP_SOURCE_DIR}/include
					 # ${libKinebenUtils_SOURCE_DIR}/include
                     #${libDecisionTree_SOURCE_DIR}/include
                     #${libMultiKeyframe_SOURCE_DIR}/include
					 ${libCC3D_SOURCE_DIR}/include
                     ${libIndexedMesh_SOURCE_DIR}/include
					 #${libRandomForest_SOURCE_DIR}/include
                     ${libGMorpher_SOURCE_DIR}/include
                     ${libPatchedCloud_SOURCE_DIR}/include
                     ${libGMorpher_PatchedMesh_SOURCE_DIR}/include
                     ${CUDA_TOOLKIT_ROOT_DIR}/include
					 ${OPENCV_INCLUDE_DIR}
					 ${Boost_INCLUDE_DIR}
					 ${EIGEN3_INCLUDE_DIR}					 
					 ${MEX_INCLUDE_DIR}		 
					 ${CERES_INCLUDE_DIR}				 
					 ${GLOG_INCLUDE_DIR})


#  first we compile all the common files to a lib because it will be used by two excutables
ADD_LIBRARY(GMorpher_EMLib STATIC ${SRC_COMMON} ${COMPILED_KERNELS})

# then we compile and link both executables
LINK_DIRECTORIES(	${CUDA_TOOLKIT_ROOT_DIR}/lib/x64 
					${OpenCV_DIR}/install/x64/vc${VC_VER}/lib
					${MEX_LIB_DIR}
					${SUITESPARSE_LIB_DIR})	 

LINK_LIBRARIES(  GMorpher
				 GMorpher_EMLib
				 GMorpher_PatchedMesh
				 #MultiKeyframe
				 #RandomForest
                 PatchedCloud
                 IndexedMesh
                 cuda cudart  
				 debug opencv_core2413d debug opencv_highgui2413d debug opencv_imgproc2413d
				 optimized opencv_core2413 optimized opencv_highgui2413 optimized opencv_imgproc2413			
                )
			
ADD_EXECUTABLE( GMorpherBone_3DEMICP    src/source_surface.cpp   ${KERNEL_SRC})	#  src/source_hough.cpp		src/source_rangeTry.cpp        src/source_FaustEx.cpp  src/source_forestCeres.cpp     src/source_CVTfromSur.cpp 

SET_TARGET_PROPERTIES(GMorpherBone_3DEMICP PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcpmt.lib /NODEFAULTLIB:libcpmtd.lib /NODEFAULTLIB:libcmtd.lib /NODEFAULTLIB:libcmt.lib")


