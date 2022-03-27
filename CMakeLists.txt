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
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

add_executable(renderer main.cpp)
add_executable(visualize_sdf visualize_sdf.cpp)

find_package(OpenCV REQUIRED)
if (OpenCV_FOUND)
	include_directories( ${OpenCV_INCLUDE_DIRS} )
	add_executable(ssim ssim.cpp)
	target_link_libraries(ssim ${OpenCV_LIBS} )
else()
	message("OpenCV not found")
endif()

target_compile_features(renderer PRIVATE cxx_std_17)
