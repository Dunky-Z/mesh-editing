
cmake_minimum_required (VERSION 3.8)

project ("human-reconstruction")

#PMP
find_package(PMP)

#Eigen
find_package(Eigen3 REQUIRED)

#SURFACE MESH
include_directories("$ENV{Surface_Mesh_DIR}/include")

include_directories("./include")
LINK_DIRECTORIES("./Lib")
#include_directories("D:\Lib\OpenGL\include")

#find_package(GLFW REQUIRED)

add_executable(human-reconstruction src/testpmp.cpp )
target_link_libraries(human-reconstruction pmp)
target_link_libraries(human-reconstruction Eigen3::Eigen)
target_link_libraries(human-reconstruction debug "$ENV{Surface_Mesh_DIR}/lib/surface_mesh_debug.lib")
#target_link_libraries(3D-HumanBodyshape debug "D:/ITabc/ITabc/BodyReconstruction/BodyReconstruction/libs/python37_d.lib")
