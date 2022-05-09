cmake_minimum_required(VERSION 3.1)
project(qp_solver_collection)


# Options for solvers
option(DEFAULT_ENABLE_ALL "Enable all QP solvers by default" OFF)
set(DEFAULT_ENABLE_VALUE ${DEFAULT_ENABLE_ALL})

option(ENABLE_QLD "Enable QLD" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_QUADPROG "Enable QuadProg" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_LSSOL "Enable LSSOL" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_JRLQP "Enable JRLQP" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_QPOASES "Enable qpOASES" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_OSQP "Enable OSQP" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_NASOQ "Enable NASOQ" ${DEFAULT_ENABLE_VALUE})

# Option to switch between Catkin and pure CMake build system
option(QP_SOLVER_COLLECTION_STANDALONE "It ON, build  without requiring Catkin and roscpp" OFF)

if(NOT QP_SOLVER_COLLECTION_STANDALONE)
  find_package(catkin REQUIRED COMPONENTS
    roscpp
    )
  add_compile_options(-std=c++17)
else()
  if(CMAKE_VERSION VERSION_LESS "3.16")
    message(FATAL_ERROR "To enable QP_SOLVER_COLLECTION_STANDALONE option you need at least CMake 3.16")
  endif()
  option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)
  option(BUILD_TESTING "Build test" OFF)
  if(BUILD_TESTING)
    enable_testing()
  endif()
  include(GNUInstallDirs)
  # Place output artifacts all in the same directory to simplify running tests on Windows
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
endif()

# To build shared libraries in Windows, we set CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS to TRUE.
# See https://cmake.org/cmake/help/v3.4/variable/CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS.html
# See https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

macro(add_qp_solver _SOLVER_NAME _PACKAGE_NAME)
  if(${ENABLE_${_SOLVER_NAME}})
    if(NOT "${_PACKAGE_NAME}" STREQUAL "")
      find_package(${_PACKAGE_NAME} REQUIRED)
    endif()
  endif()
  message(NOTICE "- ENABLE_${_SOLVER_NAME}: ${ENABLE_${_SOLVER_NAME}}")
endmacro()

# eigen-qld
add_qp_solver(QLD eigen-qld)

# eigen-quadprog
add_qp_solver(QUADPROG eigen-quadprog)

# eigen-lssol
add_qp_solver(LSSOL eigen-lssol)

# JRLQP
add_qp_solver(JRLQP jrl-qp)

# qpOASES
add_qp_solver(QPOASES "")

# OsqpEigen
add_qp_solver(OSQP OsqpEigen)

# NASOQ
add_qp_solver(NASOQ nasoq)

if(NOT QP_SOLVER_COLLECTION_STANDALONE)
  catkin_package(
    CATKIN_DEPENDS
    roscpp
    DEPENDS EIGEN3
    INCLUDE_DIRS include
    LIBRARIES QpSolverCollection
    CFG_EXTRAS qp_solver_options.cmake
    )
endif()

set(QP_SOLVER_OPTIONS_HEADER_FILE "${PROJECT_SOURCE_DIR}/include/qp_solver_collection/QpSolverOptions.h")
configure_file("${QP_SOLVER_OPTIONS_HEADER_FILE}.in" "${QP_SOLVER_OPTIONS_HEADER_FILE}")

set(QP_SOLVER_OPTIONS_CMAKE_FILE "${PROJECT_SOURCE_DIR}/cmake/qp_solver_options.cmake")
configure_file("${QP_SOLVER_OPTIONS_CMAKE_FILE}.in" "${QP_SOLVER_OPTIONS_CMAKE_FILE}")

add_subdirectory(src)

if(QP_SOLVER_COLLECTION_STANDALONE)
  install(EXPORT ${PROJECT_NAME}
    FILE ${PROJECT_NAME}Targets.cmake
    DESTINATION lib/cmake/${PROJECT_NAME}
    NAMESPACE ${PROJECT_NAME}::
    )

  # Write initial config cmake file
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in"
"
\@PACKAGE_INIT\@

include(\"\${CMAKE_CURRENT_LIST_DIR}/\@PROJECT_NAME\@Targets.cmake\")
include(\"\${CMAKE_CURRENT_LIST_DIR}/qp_solver_options.cmake\")

set(\@PROJECT_NAME\@_LIBRARIES \@PROJECT_NAME\@::QpSolverCollection)
")

  include(CMakePackageConfigHelpers)
  # generate the config file that is includes the exports
  configure_package_config_file(${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

  # generate the version file for the config file
  # Extract version numbers from package.xml
  file(READ package.xml PACKAGE_XML)
  string(REGEX MATCH "<version>[0-9]+\\.[0-9]+\\.[0-9]+</version>" DIRTY_VERSION_STRING ${PACKAGE_XML})
  string(REGEX REPLACE "^<version>([0-9]+)\\.([0-9]+)\\.([0-9]+)</version>$" "\\1" OSC_MAJOR_VERSION "${DIRTY_VERSION_STRING}")
  string(REGEX REPLACE "^<version>([0-9]+)\\.([0-9]+)\\.([0-9]+)</version>$" "\\2" QSC_MINOR_VERSION "${DIRTY_VERSION_STRING}")
  string(REGEX REPLACE "^<version>([0-9]+)\\.([0-9]+)\\.([0-9]+)</version>$" "\\3" QSC_PATCH_VERSION "${DIRTY_VERSION_STRING}")
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${QSC_MAJOR_VERSION}.${QSC_MINOR_VERSION}.${QSC_PATCH_VERSION}"
    COMPATIBILITY AnyNewerVersion
    )

  # install config files
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    "${QP_SOLVER_OPTIONS_CMAKE_FILE}"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    )
endif()

if(CATKIN_ENABLE_TESTING OR BUILD_TESTING)
  add_subdirectory(tests)
endif()

OPTION(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)
if(INSTALL_DOCUMENTATION)
  add_subdirectory(doc)
endif()
