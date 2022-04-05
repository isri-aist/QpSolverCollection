cmake_minimum_required(VERSION 3.1)
project(qp_solver_collection)

add_compile_options(-std=c++17)

# Options for solvers
option(DEFAULT_ENABLE_ALL "Enable all QP solvers by default" OFF)
set(DEFAULT_ENABLE_VALUE ${DEFAULT_ENABLE_ALL})

option(ENABLE_QLD "Enable QLD" ON)
option(ENABLE_QUADPROG "Enable QuadProg" ON)
option(ENABLE_LSSOL "Enable LSSOL" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_JRLQP "Enable JRLQP" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_QPOASES "Enable qpOASES" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_OSQP "Enable OSQP" ${DEFAULT_ENABLE_VALUE})
option(ENABLE_NASOQ "Enable NASOQ" ${DEFAULT_ENABLE_VALUE})

find_package(catkin REQUIRED COMPONENTS
  message_generation
  roscpp
  )

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

macro(add_qp_solver _SOLVER_NAME _PACKAGE_NAME)
  if(${ENABLE_${_SOLVER_NAME}})
    if(NOT "${_PACKAGE_NAME}" STREQUAL "")
      find_package(${_PACKAGE_NAME} REQUIRED)
    endif()
    add_definitions(-DENABLE_${_SOLVER_NAME}=1)
  else()
    add_definitions(-DENABLE_${_SOLVER_NAME}=0)
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

catkin_package(
  CATKIN_DEPENDS
  roscpp
  DEPENDS EIGEN3
  INCLUDE_DIRS include
  LIBRARIES QpSolverCollection
  )

include_directories(include ${catkin_INCLUDE_DIRS})

add_subdirectory(src)

if(CATKIN_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

OPTION(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)
if(INSTALL_DOCUMENTATION)
  add_subdirectory(doc)
endif()
