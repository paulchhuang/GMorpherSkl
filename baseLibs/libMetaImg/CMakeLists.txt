project (libMetaImg)

cmake_minimum_required(VERSION 2.8)
#-------------------- OpenMP part-------------------------
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
set( METAIMG_SRC src/MetaImg.cpp)

find_package(OpenCV REQUIRED)

INCLUDE_DIRECTORIES( include 
						${libMVLib_SOURCE_DIR}/include
						${OpenCV_INCLUDE_DIRS} 
						${EIGEN3_INCLUDE_DIR})

if(MSVC_VERSION STREQUAL "1800")
	set(VC_VER "12")
elseif(MSVC_VERSION STREQUAL "1900")
	set(VC_VER "14")
elseif(MSVC_VERSION STREQUAL "1910")
	set(VC_VER "15")
endif()
LINK_DIRECTORIES(${OpenCV_DIR}/install/x64/vc${VC_VER}/lib)

LINK_LIBRARIES(	MVLib 
				debug opencv_core2413d debug opencv_highgui2413d debug opencv_imgproc2413d  
				optimized opencv_core2413 optimized opencv_highgui2413 optimized opencv_imgproc2413 )


#ADD_DEFINITIONS(-DINDEXEDMESH_EXPORTS -D_CRT_NOFORCE_MANIFEST)
ADD_LIBRARY( MetaImg STATIC ${METAIMG_SRC} )


# ##############################
# Platform Specific
IF (WIN32)
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")
ENDIF(WIN32)


if ( UNIX )
	SET_TARGET_PROPERTIES( MetaImg PROPERTIES COMPILE_FLAGS -fPIC)
endif( UNIX )




