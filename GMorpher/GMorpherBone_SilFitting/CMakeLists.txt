project (GMorpherBone_SIL)

cmake_minimum_required(VERSION 2.8)

#-------------------- OpenMP part-------------------------
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)	
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

FIND_PACKAGE(Boost COMPONENTS system filesystem program_options REQUIRED)
if(Boost_FOUND)	
	MESSAGE("Boost_INCLUDE: " ${Boost_INCLUDE_DIR})
	MESSAGE("Boost_LINK: " ${Boost_LIBRARY_DIRS} )
endif()

find_package(OpenCV REQUIRED)

# -------------------------

if(MSVC_VERSION STREQUAL "1800")
	set(VC_VER "12")
elseif(MSVC_VERSION STREQUAL "1900")
	set(VC_VER "14")
elseif(MSVC_VERSION STREQUAL "1910")
	set(VC_VER "15")
endif()
# -------------------------
# THE SOURCES

# the mesh/mesh app
set( SRC
     src/EnergyTerm_SilFit.cpp
     src/source.cpp )


INCLUDE_DIRECTORIES( include
                     ${libCC3D_SOURCE_DIR}/include
                     ${libIndexedMesh_SOURCE_DIR}/include
                     ${libMetaImg_SOURCE_DIR}/include
                     ${libGMorpher_SOURCE_DIR}/include
                     ${libPatchedCloud_SOURCE_DIR}/include
                     ${libGMorpher_PatchedMesh_SOURCE_DIR}/include
                     ${libMVLib_SOURCE_DIR}/include
                     ${libRendererOffScreen_SOURCE_DIR}/include
					 ${OPENCV_INCLUDE_DIRS}
					 ${Boost_INCLUDE_DIR}
					 ${OPENGL_INCLUDE_DIR}
					 ${CERES_INCLUDE_DIR}	
					 ${GLOG_INCLUDE_DIR}
					 ${EIGEN3_INCLUDE_DIR}
					 ${MEX_INCLUDE_DIR})

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")

LINK_DIRECTORIES(${OpenCV_DIR}/install/x64/vc${VC_VER}/lib
				 ${OPENGL_LIB_DIR}	
				 ${MEX_LIB_DIR}
				 ${SUITESPARSE_LIB_DIR})

LINK_LIBRARIES( OpenGL32 GlU32 glew32 freeglut
				GMorpher 
                IndexedMesh 
                MVLib 
                MetaImg 
                RendererOffScreen 
                GMorpher_PatchedMesh 
                PatchedCloud
				debug opencv_core2413d debug opencv_highgui2413d debug opencv_imgproc2413d
				optimized opencv_core2413 optimized opencv_highgui2413 optimized opencv_imgproc2413)

ADD_EXECUTABLE( GMorpherBone_SIL ${SRC} )
#SET_TARGET_PROPERTIES(GMorpherBone_SIL PROPERTIES LINK_FLAGS "/NODEFAULTLIB:MSVCRTD.lib /NODEFAULTLIB:MSVCRT.lib /NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:msvcprtd.lib")
SET_TARGET_PROPERTIES(GMorpherBone_ImgBlob PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcpmt.lib /NODEFAULTLIB:libcpmtd.lib /NODEFAULTLIB:libcmtd.lib /NODEFAULTLIB:libcmt.lib")

#SET_TARGET_PROPERTIES(GMorpherBone_SIL PROPERTIES LINK_FLAGS "/NODEFAULTLIB:MSVCRTD.lib /NODEFAULTLIB:MSVCRT.lib /NODEFAULTLIB:libcmtd.lib /NODEFAULTLIB:msvcprtd.lib")