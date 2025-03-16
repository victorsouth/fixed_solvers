# header only library
# for install, dependencies and correct include path
project(fixed_solvers)

set(CMAKE_CXX_STANDARD 17)#23?
set(CMAKE_CXX_STANDARD_REQUIRED True)



find_package(Eigen3)

set(HEADERS_LS
    fixed/line_search/divider.h
    fixed/line_search/golden_section.h
)
set(HEADERS
    fixed/array_ext.h
    fixed/fixed.h
    fixed/fixed_bisection.h
    fixed/fixed_linear_solver.h
    fixed/fixed_nonlinear_solver.h
    fixed/fixed_optimizer.h
    fixed/fixed_system.h
)

if("3.19.0" VERSION_LESS ${CMAKE_VERSION})
    add_library(${PROJECT_NAME} INTERFACE ${HEADERS} ${HEADERS_LS})
else()
    add_library(${PROJECT_NAME} INTERFACE)
endif()
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_include_directories(${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/fixed_solvers>
)
add_library(fixed_solvers::fixed_solvers ALIAS fixed_solvers)

option(FIXED_SOLVERS_INSTALL "" ON)

if(FIXED_SOLVERS_INSTALL)
    include(CMakePackageConfigHelpers)
    include(GNUInstallDirs)
    install(TARGETS fixed_solvers
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        EXPORT fixed_solversTargets
    )
install(EXPORT fixed_solversTargets NAMESPACE fixed_solvers::
    FILE fixed_solversTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fixed_solvers
)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/fixed_solversConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fixed_solvers"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/fixed_solversConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fixed_solvers
)

install(FILES ${HEADERS} DESTINATION include/fixed_solvers/fixed)
install(FILES ${HEADERS_LS} DESTINATION include/fixed_solvers/fixed/line_search)
endif()

option(FIXED_SOLVERS_BUILD_TESTS "" OF)

if(FIXED_SOLVERS_BUILD_TESTS)

find_package(Threads)
find_package(GTest REQUIRED)

add_executable(fs_tests testing/test_main.cpp testing/bisection.cpp testing/optimizer.cpp)
target_link_libraries(fs_tests fixed_solvers::fixed_solvers GTest::gtest)
install(TARGETS fs_tests DESTINATION bin)


endif()
