/* Author: Masaki Murooka */

#pragma once

#include <chrono>
#include <fstream>
#include <memory>

#include <Eigen/SparseCore>

#include <qp_solver_collection/QpSolverOptions.h>

#include <ros/console.h>

namespace Eigen
{
class QLDDirect;
class QuadProgDense;
class LSSOL_QP;
} // namespace Eigen

namespace jrl
{
namespace qp
{
class GoldfarbIdnaniSolver;
} // namespace qp
} // namespace jrl

namespace qpOASES
{
class SQProblem;
} // namespace qpOASES

namespace OsqpEigen
{
class Solver;
} // namespace OsqpEigen

namespace QpSolverCollection
{
/** \brief QP solver type. */
enum class QpSolverType
{
  Uninitialized = -1,
  QLD = 0,
  QuadProg,
  LSSOL,
  JRLQP,
  qpOASES,
  OSQP,
  NASOQ
};

/*! \brief Convert std::string to QpSolverType. */
QpSolverType strToQpSolverType(const std::string & qp_solver_type);
} // namespace QpSolverCollection

namespace std
{
using QpSolverType = QpSolverCollection::QpSolverType;

inline string to_string(QpSolverType qp_solver_type)
{
  switch(qp_solver_type)
  {
    case QpSolverType::QLD:
      return "QpSolverType::QLD";
    case QpSolverType::QuadProg:
      return "QpSolverType::QuadProg";
    case QpSolverType::LSSOL:
      return "QpSolverType::LSSOL";
    case QpSolverType::JRLQP:
      return "QpSolverType::JRLQP";
    case QpSolverType::qpOASES:
      return "QpSolverType::qpOASES";
    case QpSolverType::OSQP:
      return "QpSolverType::OSQP";
    case QpSolverType::NASOQ:
      return "QpSolverType::NASOQ";
    default:
      ROS_ERROR_STREAM("[QpSolverType] Unsupported value: " << std::to_string(static_cast<int>(qp_solver_type)));
  }

  return "";
}
} // namespace std

namespace QpSolverCollection
{
/** \brief Class of QP coefficient.

    \todo Support both-sided inequality constraints (i.e., \f$\boldsymbol{d}_{lower} \leq \boldsymbol{C} \boldsymbol{x}
   \leq \boldsymbol{d}_{upper}\f$). QLD, QuadProg, and NASOQ support only one-sided constraints, while LSSOL, JRLQP,
   QPOASES, and OSQP support both-sided constraints.
*/
class QpCoeff
{
public:
  /** \brief Constructor. */
  QpCoeff() {}

  /** \brief Setup the coefficients with filling zero.
      \param dim_var dimension of decision variable
      \param dim_eq dimension of equality constraint
      \param dim_ineq dimension of inequality constraint
  */
  void setup(int dim_var, int dim_eq, int dim_ineq);

  /** \brief Print information. */
  void printInfo(bool verbose = false, const std::string & header = "") const;

  /** \brief Dump coefficients. */
  void dump(std::ofstream & ofs) const;

public:
  /** \brief Dimensions.
      @{
  */
  int dim_var_ = 0;
  int dim_eq_ = 0;
  int dim_ineq_ = 0;
  /** @} */

  /** \brief Coefficients.
      @{
  */
  Eigen::MatrixXd obj_mat_;
  Eigen::VectorXd obj_vec_;
  Eigen::MatrixXd eq_mat_;
  Eigen::VectorXd eq_vec_;
  Eigen::MatrixXd ineq_mat_;
  Eigen::VectorXd ineq_vec_;
  Eigen::VectorXd x_min_;
  Eigen::VectorXd x_max_;
  /** @} */
};

/** \brief Virtual class of QP solver. */
class QpSolver
{
public:
  using clock = typename std::conditional<std::chrono::high_resolution_clock::is_steady,
                                          std::chrono::high_resolution_clock,
                                          std::chrono::steady_clock>::type;

public:
  /** \brief Constructor. */
  QpSolver() {}

  /** \brief Print information. */
  void printInfo(bool verbose = false, const std::string & header = "") const;

  /** \brief Solve QP.
      \param dim_var dimension of decision variable
      \param dim_eq dimension of equality constraint
      \param dim_ineq dimension of inequality constraint
      \param Q objective matrix (LSSOL requires non-const for Q)
      \param c objective vector
      \param A equality constraint matrix
      \param b equality constraint vector
      \param C inequality constraint matrix
      \param d inequality constraint vector
      \param x_min lower bound
      \param x_max upper bound

      QP is formulated as follows:
      \f{align*}{
      & min_{\boldsymbol{x}} \ \frac{1}{2}{\boldsymbol{x}^T \boldsymbol{Q} \boldsymbol{x}} + {\boldsymbol{c}^T
     \boldsymbol{x}} \\
      & s.t. \ \ \boldsymbol{A} \boldsymbol{x} = \boldsymbol{b} \nonumber \\
      & \phantom{s.t.} \ \ \boldsymbol{C} \boldsymbol{x} \leq \boldsymbol{d} \nonumber \\
      & \phantom{s.t.} \ \ \boldsymbol{x}_{min} \leq \boldsymbol{x} \leq \boldsymbol{x}_{max} \nonumber
      \f}
  */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) = 0;

  /** \brief Solve QP.
      \param qp_coeff QP coefficient
  */
  virtual Eigen::VectorXd solve(QpCoeff & qp_coeff);

public:
  /** \brief QP solver type. */
  QpSolverType type_ = QpSolverType::Uninitialized;

  /** \brief Whether it failed to solve the QP. */
  bool solve_failed_ = false;
};

#if ENABLE_QLD
/** \brief QP solver QLD. */
class QpSolverQld : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverQld();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

protected:
  std::unique_ptr<Eigen::QLDDirect> qld_;
};
#endif

#if ENABLE_QUADPROG
/** \brief QP solver QuadProg. */
class QpSolverQuadprog : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverQuadprog();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

protected:
  std::unique_ptr<Eigen::QuadProgDense> quadprog_;
};
#endif

#if ENABLE_LSSOL
/** \brief QP solver LSSOL. */
class QpSolverLssol : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverLssol();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

protected:
  std::unique_ptr<Eigen::LSSOL_QP> lssol_;
};
#endif

#if ENABLE_JRLQP
/** \brief QP solver JRLQP. */
class QpSolverJrlqp : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverJrlqp();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

protected:
  std::unique_ptr<jrl::qp::GoldfarbIdnaniSolver> jrlqp_;
};
#endif

#if ENABLE_QPOASES
/** \brief QP solver qpOASES.
    \todo Support an efficient interface (QProblemB) dedicated to QP with only box constraints.
*/
class QpSolverQpoases : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverQpoases();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

public:
  int n_wsr_ = 10000;

  /** \brief Whether to initialize each time instead of doing a warm start.

      \note Warm start did not give good results.
  */
  bool force_initialize_ = true;

protected:
  std::unique_ptr<qpOASES::SQProblem> qpoases_;
};
#endif

#if ENABLE_OSQP
/** \brief QP solver OSQP.
    \todo Set without going through a dense matrix.
 */
class QpSolverOsqp : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverOsqp();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

public:
  /** \brief Whether to initialize each time instead of doing a warm start.

      \note Warm start did not give good results.
  */
  bool force_initialize_ = true;

protected:
  std::unique_ptr<OsqpEigen::Solver> osqp_;

  Eigen::SparseMatrix<double> Q_sparse_;
  Eigen::VectorXd c_;
  Eigen::SparseMatrix<double> AC_with_bound_sparse_;
  Eigen::VectorXd bd_with_bound_min_;
  Eigen::VectorXd bd_with_bound_max_;

  double sparse_duration_ = 0; // [ms]
};
#endif

#if ENABLE_NASOQ
/** \brief QP solver NASOQ.
    \todo Set without going through a dense matrix.
*/
class QpSolverNasoq : public QpSolver
{
public:
  /** \brief Constructor. */
  QpSolverNasoq();

  /** \brief Solve QP. */
  virtual Eigen::VectorXd solve(int dim_var,
                                int dim_eq,
                                int dim_ineq,
                                Eigen::Ref<Eigen::MatrixXd> Q,
                                const Eigen::Ref<const Eigen::VectorXd> & c,
                                const Eigen::Ref<const Eigen::MatrixXd> & A,
                                const Eigen::Ref<const Eigen::VectorXd> & b,
                                const Eigen::Ref<const Eigen::MatrixXd> & C,
                                const Eigen::Ref<const Eigen::VectorXd> & d,
                                const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                const Eigen::Ref<const Eigen::VectorXd> & x_max) override;

protected:
  Eigen::SparseMatrix<double, Eigen::ColMajor, int> Q_sparse_;
  Eigen::SparseMatrix<double, Eigen::ColMajor, int> A_sparse_;
  Eigen::SparseMatrix<double, Eigen::ColMajor, int> C_with_bound_sparse_;

  double sparse_duration_ = 0; // [ms]
};
#endif

/** \brief Check whether QP solver is enabled.
    \param qp_solver_type QP solver type
 */
bool isQpSolverEnabled(const QpSolverType & qp_solver_type);

/** \brief Allocate the specified QP solver.
    \param qp_solver_type QP solver type
 */
std::shared_ptr<QpSolver> allocateQpSolver(const QpSolverType & qp_solver_type);
} // namespace QpSolverCollection
