project(OpenFBX LANGUAGES CXX C)

#find_package(miniz REQUIRED)

add_library(${PROJECT_NAME} ofbx.cpp libdeflate.c)
#target_link_libraries(${PROJECT_NAME} PRIVATE miniz::miniz)

target_include_directories(${PROJECT_NAME}
        PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_compile_definitions(${PROJECT_NAME} PRIVATE _LARGEFILE64_SOURCE)
set_target_properties(${PROJECT_NAME} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
