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

TEST(TestQpSolversEnabled, QLD)
{
  checkOneSolver(QpSolverType::QLD);
}

TEST(TestQpSolversEnabled, QuadProg)
{
  checkOneSolver(QpSolverType::QuadProg);
}

TEST(TestQpSolversEnabled, LSSOL)
{
  checkOneSolver(QpSolverType::LSSOL);
}

TEST(TestQpSolversEnabled, JRLQP)
{
  checkOneSolver(QpSolverType::JRLQP);
}

TEST(TestQpSolversEnabled, qpOASES)
{
  checkOneSolver(QpSolverType::qpOASES);
}

TEST(TestQpSolversEnabled, OSQP)
{
  checkOneSolver(QpSolverType::OSQP);
}

TEST(TestQpSolversEnabled, NASOQ)
{
  checkOneSolver(QpSolverType::NASOQ);
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
