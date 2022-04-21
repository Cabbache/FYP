cmake_minimum_required(VERSION 3.16.3)
project(
	Renderer
	VERSION
		1.0
	LANGUAGES
		CXX
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(SDL2 REQUIRED)
find_package(OpenCV REQUIRED)

add_executable(visualize_sdf visualize_sdf.cpp)

if (OpenCV_FOUND)
	include_directories( ${OpenCV_INCLUDE_DIRS} )
	add_executable(ssim ssim.cpp)
	add_executable(psnr psnr.cpp)
	target_link_libraries(ssim ${OpenCV_LIBS} )
	target_link_libraries(psnr ${OpenCV_LIBS} )
else()
	message("OpenCV not found")
endif()

add_executable(renderer main.cpp headers/KDTreeCPU.cpp headers/mesh.cpp headers/KDTreeNode.cpp)
target_compile_features(renderer PRIVATE cxx_std_17)
target_link_libraries(renderer ${SDL2_LIBRARIES})

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
