
cmake_minimum_required (VERSION 3.8)

project ("mesh-editing")

#PMP
find_package(PMP)

#Eigen
find_package(Eigen3 REQUIRED)

#SURFACE MESH
include_directories("$ENV{Surface_Mesh_DIR}/include")
message($ENV{Surface_Mesh_DIR})


#find_package(GLFW REQUIRED)

add_executable(mesh-editing src/testpmp.cpp )
target_link_libraries(mesh-editing pmp)
target_link_libraries(mesh-editing Eigen3::Eigen)
#target_link_libraries(mesh-editing debug "$ENV{Surface_Mesh_DIR}/lib/SurfaceMesh.lib")
#target_link_libraries(mesh-editing debug "D:/ITabc/ITabc/BodyReconstruction/BodyReconstruction/libs/python37_d.lib")
