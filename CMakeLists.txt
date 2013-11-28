cmake_minimum_required (VERSION 2.6) 
project (ContourSpectrum) 

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)
#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)
#########################################################
# FIND GLEW
#########################################################
find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARY_DIRS})
add_definitions(${GLEW_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)
#########################################################

if(NOT VTK_BINARY_DIR)
  find_package(VTK 
 # COMPONENTS
 #   vtkCommonCore
 #   vtkCommonDataModel
 #   vtkIOLegacy
  )
  include(${VTK_USE_FILE})
endif()



add_executable (ContourSpectrum 
				ContourSpectrum.cpp 
				DataMgr.cpp 
				Renderer.cpp
				
				ContourSpectrum.h
				DataMgr.h
				Renderer.h
				my_vector_types.h
				my_vector_functions.h
				) 
target_link_libraries(ContourSpectrum ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${VTK_LIBRARIES} ${GLEW_LIBRARY})

