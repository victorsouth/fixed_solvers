# header only library
# for install and correct include path
project(fixed_solvers)

set(HEADERS_LS
fixed/line_search/divider.h
fixed/line_search/golden_section.h
)
set(HEADERS
fixed/array_ext.h
fixed/fixed.h
fixed/fixed_linear_solver.h
fixed/fixed_nonlinear_solver.h
fixed/fixed_system.h
)
include(CMakePackageConfigHelpers)

if("3.19.0" VERSION_LESS ${CMAKE_VERSION})
add_library(${PROJECT_NAME} INTERFACE ${HEADERS} ${HEADERS_LS})
else()
endif()
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME}
                             INTERFACE
                        $<INSTALL_INTERFACE:include/fixed_solvers>
)
add_library(fixed_solvers::fixed_solvers ALIAS fixed_solvers)
install(TARGETS fixed_solvers
        DESTINATION lib
        EXPORT fixed_solversTargets
)
install(EXPORT fixed_solversTargets NAMESPACE fixed_solvers::
  FILE fixed_solversTargets.cmake
  DESTINATION lib/cmake/fixed_solvers
)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/fixed_solversConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/fixed_solvers"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/fixed_solversConfig.cmake
  DESTINATION lib/cmake/fixed_solvers
  )

install(FILES ${HEADERS} DESTINATION include/fixed_solvers/fixed)
install(FILES ${HEADERS_LS} DESTINATION include/fixed_solvers/fixed/line_search)
