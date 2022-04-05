set(QpSolverCollection_gtest_list
  TestQpSolversEnabled
  TestSampleQP
  )

foreach(NAME IN LISTS QpSolverCollection_gtest_list)
  catkin_add_gtest(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} QpSolverCollection)
endforeach()
