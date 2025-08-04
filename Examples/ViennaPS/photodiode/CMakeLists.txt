find_package(OpenMP)
find_package(ViennaPS 3.5.0)

if(NOT TARGET ViennaTools::ViennaPS)
  message(STATUS "ViennaPS not found, skipping example")
  return()
endif()


project(photodiode LANGUAGES CXX)

add_executable("${PROJECT_NAME}" "${PROJECT_NAME}.cpp")
target_link_libraries("${PROJECT_NAME}" PRIVATE FiDiTi ViennaTools::ViennaPS)

configure_file(config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)

add_dependencies(FiDiTi_Examples ${PROJECT_NAME})
