cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(LBC CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Eigen3
find_package(Eigen3 REQUIRED)

find_package(OpenMP)

set(${PROJECT_NAME}_SOURCES
	LBCSolver.h
    LBCSolver.cpp
    DataSetup.h
    DataSetup.cpp
    common.h
)

add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})

target_link_libraries(${PROJECT_NAME}
    PUBLIC
    Eigen3::Eigen
)

if(OpenMP_FOUND)
	target_compile_options(
		${PROJECT_NAME}
		PUBLIC
		$<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>
	)
	target_link_libraries(
	${PROJECT_NAME}
	PUBLIC
	${OpenMP_CXX_LIBRARIES}
	)
endif()

target_include_directories(
	${PROJECT_NAME}
	INTERFACE
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)