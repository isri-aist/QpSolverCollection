/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <qp_solver_collection/QpSolverCollection.h>

using QpSolverCollection::QpSolverType;

void checkOneSolver(QpSolverType qp_solver_type)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(qp_solver_type));
  EXPECT_TRUE(QpSolverCollection::allocateQpSolver(qp_solver_type));
}

TEST(TestQpSolversEnabled, Any)
{
  checkOneSolver(QpSolverType::Any);
}

#if ENABLE_QLD
TEST(TestQpSolversEnabled, QLD)
{
  checkOneSolver(QpSolverType::QLD);
}
#endif

#if ENABLE_QUADPROG
TEST(TestQpSolversEnabled, QuadProg)
{
  checkOneSolver(QpSolverType::QuadProg);
}
#endif

#if ENABLE_LSSOL
TEST(TestQpSolversEnabled, LSSOL)
{
  checkOneSolver(QpSolverType::LSSOL);
}
#endif

#if ENABLE_JRLQP
TEST(TestQpSolversEnabled, JRLQP)
{
  checkOneSolver(QpSolverType::JRLQP);
}
#endif

#if ENABLE_QPOASES
TEST(TestQpSolversEnabled, qpOASES)
{
  checkOneSolver(QpSolverType::qpOASES);
}
#endif

#if ENABLE_OSQP
TEST(TestQpSolversEnabled, OSQP)
{
  checkOneSolver(QpSolverType::OSQP);
}
#endif

#if ENABLE_NASOQ
TEST(TestQpSolversEnabled, NASOQ)
{
  checkOneSolver(QpSolverType::NASOQ);
}
#endif

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
