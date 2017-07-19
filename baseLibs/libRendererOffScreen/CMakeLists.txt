project (libRendererOffScreen)

cmake_minimum_required(VERSION 2.8)


set( SRC src/GLerrors.cpp	#src/OpenGLContextMS.cpp
         src/Renderer.cpp)		 

INCLUDE_DIRECTORIES( include ${OPENGL_INCLUDE_DIR})

LINK_DIRECTORIES(${OPENGL_LIB_DIR})
LINK_LIBRARIES(OpenGL32 GlU32 glew32 freeglut)
#SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3 -march=athlon64 -m3dnow -mmmx -msse -msse2 -mfpmath=sse  -Wextra")
#SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O0 -g -D GL_GLEXT_PROTOTYPES")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")

ADD_LIBRARY( RendererOffScreen STATIC ${SRC} )




