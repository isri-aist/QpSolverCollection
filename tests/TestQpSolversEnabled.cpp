/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <qp_solver_collection/QpSolverCollection.h>

using QpSolverCollection::QpSolverType;

TEST(TestQpSolversEnabled, QLD)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::QLD));
}

TEST(TestQpSolversEnabled, QuadProg)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::QuadProg));
}

TEST(TestQpSolversEnabled, LSSOL)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::LSSOL));
}

TEST(TestQpSolversEnabled, JRLQP)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::JRLQP));
}

TEST(TestQpSolversEnabled, qpOASES)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::qpOASES));
}

TEST(TestQpSolversEnabled, OSQP)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::OSQP));
}

TEST(TestQpSolversEnabled, NASOQ)
{
  EXPECT_TRUE(QpSolverCollection::isQpSolverEnabled(QpSolverType::NASOQ));
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
