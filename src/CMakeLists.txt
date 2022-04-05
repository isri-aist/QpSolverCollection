add_library(QpSolverCollection
  QpSolver.cpp
  )

target_link_libraries(QpSolverCollection PUBLIC
  ${catkin_LIBRARIES})
if(${ENABLE_QLD})
  target_link_libraries(QpSolverCollection PUBLIC
    eigen-qld::eigen-qld)
endif()
if(${ENABLE_QUADPROG})
  target_link_libraries(QpSolverCollection PUBLIC
    eigen-quadprog::eigen-quadprog)
endif()
if(${ENABLE_LSSOL})
  target_link_libraries(QpSolverCollection PUBLIC
    eigen-lssol::eigen-lssol)
endif()
if(${ENABLE_JRLQP})
  target_link_libraries(QpSolverCollection PUBLIC
    jrl-qp::jrl-qp)
endif()
if(${ENABLE_QPOASES})
  target_include_directories(QpSolverCollection PUBLIC
    ${QPOASES_INCLUDE_DIR})
  target_link_directories(QpSolverCollection PUBLIC
    ${QPOASES_LIBRARY_DIR})
  target_link_libraries(QpSolverCollection PUBLIC
    qpOASES)
endif()
if(${ENABLE_OSQP})
  target_link_libraries(QpSolverCollection PUBLIC
    OsqpEigen::OsqpEigen)
endif()
if(${ENABLE_NASOQ})
  target_link_libraries(QpSolverCollection PUBLIC
    nasoq)
endif()
