project (libGMorpher)

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
set( GMORPHER_SRC
     src/RigidT.cpp
     src/GraphUtils.cpp
	 src/EnergyTerm_3DConstraint_ceres.cpp		
	 src/IEnergyTerm_ceres.cpp			
	 src/Solver_ceres.cpp
	 src/JointAnglePrior.cpp
	 #src/ISparseSolver.cpp		#src/SparseSolver_TAUCS.cpp     
	 #src/IEnergyTerm.cpp  src/Solver.cpp
	 #src/EnergyTerm_3DConstraint.cpp src/EnergyTerm_Rigidity_PAULY.cpp
	 #src/BigMatrix.cpp   
)


# -------------------------
# THE LIBRARY STUFF COMPILED WITH G++
	# compile
INCLUDE_DIRECTORIES( include 
                     ${libCC3D_SOURCE_DIR}/include
                     ${TAUCS_INC_DIR}
					 ${Boost_INCLUDE_DIR} 	
					 ${CERES_INCLUDE_DIR}				 
					 ${GLOG_INCLUDE_DIR}
					 ${EIGEN3_INCLUDE_DIR}
					 ${MEX_INCLUDE_DIR} )

# link
LINK_DIRECTORIES(${MEX_LIB_DIR} ${SUITESPARSE_LIB_DIR}) 	#${TAUCS_BIN_DIR} 

LINK_LIBRARIES(	debug ${CERES_LIB_DIR}/Debug/ceres-debug.lib 	optimized ${CERES_LIB_DIR}/Release/ceres.lib 
				debug ${GLOG_LIB_DIR}/Debug/libglog.lib 		optimized ${GLOG_LIB_DIR}/Release/libglog.lib 
				libblas liblapack libamd libcamd libccolamd libcholmod libcolamd libspqr suitesparseconfig metis)
# output vcf2cCG clapackCG libmetisCG blasCG libmat.lib libmx.lib 
ADD_LIBRARY( GMorpher STATIC ${GMORPHER_SRC} )