project(green CXX CUDA Fortran)
# CUDA runtime
find_package(CUDAToolkit REQUIRED)
#CGAL
find_package(CGAL REQUIRED)
#SPDLog
find_package(spdlog REQUIRED)
# LAPACK
find_package(LAPACK REQUIRED)
add_library(cuda_green cuda_impl.cu)
set_target_properties(cuda_green PROPERTIES CUDA_ARCHITECTURES native)
target_link_libraries(cuda_green
  CUDA::cudart
  Eigen3::Eigen
  CGAL::CGAL
  )

file(GLOB source *.cc *.cpp *.f90)
add_library(green ${source})
target_link_libraries(green
  igl::core
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  cuda_green
  spdlog::spdlog_header_only
  )

target_include_directories(
	green
	INTERFACE
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)
