project (libMVLib)

cmake_minimum_required(VERSION 2.8)
SET( CMAKE_BUILD_TYPE Release )

set( SRC src/BuildFilename.cpp
         src/Camera.cpp
         src/CameraSetup.cpp
         src/ImgProcess.cpp
         src/ColorConvert.cpp	# src/PNG_IO.cpp        
         src/LineDrawer.cpp)
INCLUDE_DIRECTORIES(include  ${EIGEN3_INCLUDE_DIR})


ADD_LIBRARY( MVLib STATIC ${SRC} )




