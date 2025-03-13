find_package(GTest REQUIRED)
include(GoogleTest)

set(QpSolverCollection_gtest_list
  TestQpSolversEnabled
  TestSampleQP
  )

foreach(NAME IN LISTS QpSolverCollection_gtest_list)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} GTest::gtest QpSolverCollection)
  gtest_discover_tests(${NAME})
  if(SKIP_PRIVATE_SOLVER_TEST)
    target_compile_definitions(${NAME} PRIVATE SKIP_PRIVATE_SOLVER_TEST)
  endif()
  if(FORCE_ALL_SOLVER_TEST)
    target_compile_definitions(${NAME} PRIVATE FORCE_ALL_SOLVER_TEST)
  endif()
endforeach()
