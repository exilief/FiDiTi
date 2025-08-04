#add_custom_target(FiDiTi_Examples ALL)
#
#file(GLOB entries "${CMAKE_CURRENT_LIST_DIR}/*")
#
#foreach(entry IN LISTS entries)
#    if(NOT IS_DIRECTORY "${entry}")
#        continue()
#    endif()
#
#    message(STATUS "Adding ${entry}")
#    add_subdirectory(${entry})
#endforeach()



add_custom_target(FiDiTi_Examples ALL)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/ViennaPS/photodiode")
