cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(FiDiTi
  LANGUAGES CXX C
  VERSION 0.0.0)

option(FIDITI_BUILD_EXAMPLES "Build examples" OFF)


include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 17
             CXX_EXTENSIONS OFF
             CXX_STANDARD_REQUIRED ON
             WINDOWS_EXPORT_ALL_SYMBOLS ON)


target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FiDiTi-${PROJECT_VERSION}>)


if(FIDITI_BUILD_EXAMPLES)
  message(STATUS "Building Examples")
  add_subdirectory(Examples)
endif()



include(CMakePackageConfigHelpers)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FiDiTi-${PROJECT_VERSION})

#install(FILES ...
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FiDiTi-${PROJECT_VERSION})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FiDiTi-${PROJECT_VERSION}) # FILES_MATCHING PATTERN "*.hpp")

install(EXPORT ${PROJECT_NAME}Targets
        NAMESPACE FiDiTi::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FiDiTi-${PROJECT_VERSION})

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FiDiTi-${PROJECT_VERSION})

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FiDiTi-${PROJECT_VERSION})


#export(TARGETS ${PROJECT_NAME} NAMESPACE FiDiTi:: FILE ${PROJECT_NAME}Targets.cmake)
#export(PACKAGE ${PROJECT_NAME})
