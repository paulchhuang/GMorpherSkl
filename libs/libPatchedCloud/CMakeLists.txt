project (libPatchedCloud)

cmake_minimum_required(VERSION 2.8)
#-------------------- OpenMP part-------------------------
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set( SRC src/Patching.cpp
         src/Patching_func.cpp
         src/PatchedCloud.cpp
         src/cloud_build.cpp
         src/cloud_transform.cpp
         src/cloud_interp.cpp
		 src/skl_build.cpp
		 src/skl_transform.cpp
		 src/skl_util.cpp
		 src/skl_fromCloud.cpp
		 src/Skl.cpp
		 src/RiggedMesh.cpp
)

INCLUDE_DIRECTORIES( include 
					${libCC3D_SOURCE_DIR}/include
					${libIndexedMesh_SOURCE_DIR}/include 
					${EIGEN3_INCLUDE_DIR} )

LINK_LIBRARIES( CC3D IndexedMesh )
ADD_LIBRARY( PatchedCloud STATIC ${SRC} )




