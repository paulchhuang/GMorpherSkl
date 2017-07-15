cmake_minimum_required(VERSION 2.8)

SET( LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin )
SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

SET( CMAKE_BUILD_TYPE Debug )
# SET( CMAKE_BUILD_TYPE Release )

set(THIRD_PARTY_DIR "D:/ResearchCodeMaterial/ThirdPartyLibs")

MESSAGE("CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

#-------------------- check and set up all dependencies
FIND_PACKAGE(Boost COMPONENTS system filesystem program_options REQUIRED)
if(Boost_FOUND)	
	MESSAGE("Boost_INCLUDE: " ${Boost_INCLUDE_DIR})
	MESSAGE("Boost_LINK: " ${Boost_LIBRARY_DIRS} )
endif()
FIND_PACKAGE(CGAL COMPONENTS core)
if(CGAL_FOUND)	
	MESSAGE("CGAL_DIR: " ${CGAL_DIR})
endif()

SET(CERES_INCLUDE_DIR ${THIRD_PARTY_DIR}/source/ceres-solver-1.12.0/include ${THIRD_PARTY_DIR}/source/ceres-solver-1.12.0//config)
SET(CERES_LIB_DIR	${THIRD_PARTY_DIR}/build/ceres-1_12_0/lib/)
SET(SUITESPARSE_LIB_DIR	${THIRD_PARTY_DIR}/build/suitesparse-metis-for-windows/install/lib64/
						${THIRD_PARTY_DIR}/build/suitesparse-metis-for-windows/install/lib64//lapack_blas_windows)
SET(GLOG_INCLUDE_DIR ${THIRD_PARTY_DIR}/build/glog/src/windows)
SET(GLOG_LIB_DIR ${THIRD_PARTY_DIR}/build/glog/x64)
SET(CGAL_SRC_INCLUDE ${THIRD_PARTY_DIR}/source/CGAL-4.9/include)
SET(EIGEN3_INCLUDE_DIR ${THIRD_PARTY_DIR}/source/Eigen_3_2_10 )
SET(MEX_INCLUDE_DIR "C:/Program Files/MATLAB/R2016b/extern/include")
SET(MEX_LIB_DIR "C:/Program Files/MATLAB/R2016b/extern/lib/win64/microsoft")
SET(OPENGL_INCLUDE_DIR 	"C:/Program Files (x86)/Microsoft SDKs/Windows/v7.1A/Include" 	# GL.h GLU.h
						${THIRD_PARTY_DIR}/build/glew-2.0.0/include 					# glew
						${THIRD_PARTY_DIR}/build/freeglut/include)						# glut

SET(OPENGL_LIB_DIR 	"C:/Program Files (x86)/Microsoft SDKs/Windows/v7.1A/lib/x64" 		
					${THIRD_PARTY_DIR}/build/glew-2.0.0/lib/Release/x64 
					${THIRD_PARTY_DIR}/build/freeglut/lib/x64)
 
add_subdirectory(libs)
add_subdirectory(baseLibs)
add_subdirectory(GMorpher)
# add_subdirectory(ForestTrain)