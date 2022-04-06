add_library(QpSolverCollection
  QpSolverCollection.cpp
  )

target_include_directories(QpSolverCollection PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${catkin_INCLUDE_DIRS}
  )
target_link_libraries(QpSolverCollection PUBLIC
  ${catkin_LIBRARIES}
  )
if(${ENABLE_QLD})
  target_link_libraries(QpSolverCollection PRIVATE
    eigen-qld::eigen-qld)
endif()
if(${ENABLE_QUADPROG})
  target_link_libraries(QpSolverCollection PRIVATE
    eigen-quadprog::eigen-quadprog)
endif()
if(${ENABLE_LSSOL})
  target_link_libraries(QpSolverCollection PRIVATE
    eigen-lssol::eigen-lssol)
endif()
if(${ENABLE_JRLQP})
  target_link_libraries(QpSolverCollection PRIVATE
    jrl-qp::jrl-qp)
endif()
if(${ENABLE_QPOASES})
  target_include_directories(QpSolverCollection PRIVATE
    ${QPOASES_INCLUDE_DIR})
  target_link_directories(QpSolverCollection PRIVATE
    ${QPOASES_LIBRARY_DIR})
  target_link_libraries(QpSolverCollection PRIVATE
    qpOASES)
endif()
if(${ENABLE_OSQP})
  target_link_libraries(QpSolverCollection PRIVATE
    OsqpEigen::OsqpEigen)
endif()
if(${ENABLE_NASOQ})
  target_link_libraries(QpSolverCollection PRIVATE
    nasoq)
endif()
