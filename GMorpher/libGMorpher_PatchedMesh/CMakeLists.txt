project (libGMorpher_PatchedMesh)

cmake_minimum_required(VERSION 2.8)

#-------------------- OpenMP part-------------------------
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)	
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# -------------------------
# THE SOURCES
file(GLOB SRC	"src/*.cpp")

# -------------------------
# THE LIBRARY STUFF COMPILED WITH G++
INCLUDE_DIRECTORIES(include
                    ${libGMorpher_SOURCE_DIR}/include
                    ${libPatchedCloud_SOURCE_DIR}/include
                    ${libIndexedMesh_SOURCE_DIR}/include
                    ${libCC3D_SOURCE_DIR}/include 
					${EIGEN3_INCLUDE_DIR}		
					${CERES_INCLUDE_DIR}
					${GLOG_INCLUDE_DIR}
					${Boost_INCLUDE_DIR} )
	# output
ADD_LIBRARY( GMorpher_PatchedMesh STATIC ${SRC} )





