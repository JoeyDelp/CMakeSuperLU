cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

include(GNUInstallDirs)

# Install target helper
function(superlu_target name sources)
  string(TOLOWER ${name} lcname)

  add_library(${lcname} ${sources})

  install(TARGETS ${lcname}
          EXPORT ${lcname}-targets
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(EXPORT ${lcname}-targets
          FILE CMakeSuperLUConfig.cmake
          NAMESPACE superlu::
          DESTINATION share/CMakeSuperLU/)

  set_target_properties(${lcname} PROPERTIES OUTPUT_NAME ${lcname})

  # Add include directories to target
  target_include_directories(
    ${lcname}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
           $<INSTALL_INTERFACE:include>)

  # Install include directory
  install(DIRECTORY "${PROJECT_SOURCE_DIR}/include"
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/")
endfunction()