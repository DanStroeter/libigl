cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(cageDeformation3D CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WITH_SOMIGLIANA "Compile with Somigliana coordinates (Linux only)" OFF)

# Boost
find_package(Boost REQUIRED COMPONENTS program_options)

# Eigen3
find_package(Eigen3 REQUIRED)

# IGL
set(IGL_INCLUDE_DIR "libigl/include")

add_subdirectory(LBC LBCLib)
if (WITH_SOMIGLIANA)
	add_subdirectory(somigliana SomiglianaLib)
endif()
add_subdirectory(libigl libIGL)
add_subdirectory(OpenFBX OpenFBXLib)

# Compilation using MSVC
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

set(${PROJECT_NAME}_SOURCES
	main.cpp
	globals.h
	globals.cpp
	GreenCoordinates.h
	GreenCoordinates.cpp
	MaximumLikelihoodCoordinates.h
	MaximumLikelihoodCoordinates.cpp
	MaximumEntropyCoordinates.h
	MaximumEntropyCoordinates.cpp
	LoadMesh.h
	LoadMesh.cpp
	LoadFBX.h
	LoadFBX.cpp
	InfluenceMap.h
	InfluenceMap.cpp
	Parametrization.h
	Parametrization.cpp
	WeightInterpolation.h
	WeightInterpolation.cpp
)

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})

set_target_properties(${PROJECT_NAME} PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY_RELEASE	"${CMAKE_BINARY_DIR}/build/release/bin"
	RUNTIME_OUTPUT_DIRECTORY_DEBUG		"${CMAKE_BINARY_DIR}/build/debug/bin")

set_directory_properties(PROPERTIES VS_STARTUP_PROJECT ${PROJECT_NAME})

# Compile definitions for boost
target_compile_definitions(${PROJECT_NAME} PUBLIC BOOST_PROGRAM_OPTIONS_DYN_LINK)
target_compile_definitions(${PROJECT_NAME} PUBLIC BOOST_ALL_DYN_LINK)

target_link_libraries(
	${PROJECT_NAME}
	PRIVATE
	Boost::program_options
	Boost::disable_autolinking
	Eigen3::Eigen
	LBC
	OpenFBX
)

if (WITH_SOMIGLIANA)
	target_link_libraries(
		${PROJECT_NAME}
		PRIVATE
		green
	)
	target_compile_definitions(${PROJECT_NAME} PRIVATE WITH_SOMIGLIANA)
endif()

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
    ${IGL_INCLUDE_DIR}
)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
