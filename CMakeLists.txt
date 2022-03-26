cmake_minimum_required(VERSION 3.16.3)
project(
	Renderer
	VERSION
		1.0
	LANGUAGES
		CXX
)

#default build type to release
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()


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
