/* Author: Masaki Murooka */

#define NOMINMAX

#include <limits>
#include <sstream>

#include <qp_solver_collection/QpSolverCollection.h>

// clang-format off
#if ENABLE_QLD
#include <eigen-qld/QLD.h>
#endif
#if ENABLE_QUADPROG
#include <eigen-quadprog/QuadProg.h>
#endif
#if ENABLE_LSSOL
#include <eigen-lssol/LSSOL_QP.h>
#endif
#if ENABLE_JRLQP
#include <jrl-qp/GoldfarbIdnaniSolver.h>
#include <jrl-qp/utils/enumsIO.h>
#endif
#if ENABLE_QPOASES
#include <qpOASES.hpp>
#endif
#if ENABLE_OSQP
#include <OsqpEigen/OsqpEigen.h>
#define OSQP_EIGEN_DEBUG_OUTPUT
#endif
#if ENABLE_NASOQ
#include <nasoq/nasoq_eigen.h>
#endif
// clang-format on

using namespace QpSolverCollection;

QpSolverType QpSolverCollection::strToQpSolverType(const std::string & qp_solver_type)
{
  if(qp_solver_type == "Uninitialized")
  {
    return QpSolverType::Uninitialized;
  }
  else if(qp_solver_type == "QLD")
  {
    return QpSolverType::QLD;
  }
  else if(qp_solver_type == "QuadProg")
  {
    return QpSolverType::QuadProg;
  }
  else if(qp_solver_type == "LSSOL")
  {
    return QpSolverType::LSSOL;
  }
  else if(qp_solver_type == "JRLQP")
  {
    return QpSolverType::JRLQP;
  }
  else if(qp_solver_type == "qpOASES")
  {
    return QpSolverType::qpOASES;
  }
  else if(qp_solver_type == "OSQP")
  {
    return QpSolverType::OSQP;
  }
  else if(qp_solver_type == "NASOQ")
  {
    return QpSolverType::NASOQ;
  }
  else
  {
    throw std::runtime_error("[strToQpSolverType] Unsupported QpSolverType name: " + qp_solver_type);
  }
}

void QpCoeff::setup(int dim_var, int dim_eq, int dim_ineq)
{
  dim_var_ = dim_var;
  dim_eq_ = dim_eq;
  dim_ineq_ = dim_ineq;

  obj_mat_.setZero(dim_var, dim_var);
  obj_vec_.setZero(dim_var);
  eq_mat_.setZero(dim_eq, dim_var);
  eq_vec_.setZero(dim_eq);
  ineq_mat_.setZero(dim_ineq, dim_var);
  ineq_vec_.setZero(dim_ineq);
  x_min_.setConstant(dim_var, std::numeric_limits<double>::lowest());
  x_max_.setConstant(dim_var, std::numeric_limits<double>::max());
}

void QpCoeff::printInfo(bool, const std::string & header) const
{
  QSC_INFO_STREAM(header << "dim_var: " << dim_var_ << ", dim_eq: " << dim_eq_ << ", dim_ineq: " << dim_ineq_);
}

void QpCoeff::dump(std::ofstream & ofs) const
{
  ofs << "dim_var: " << dim_var_ << std::endl;
  ofs << "dim_eq: " << dim_eq_ << std::endl;
  ofs << "dim_ineq: " << dim_ineq_ << std::endl;
  ofs << "obj_mat:\n" << obj_mat_ << std::endl;
  ofs << "obj_vec:\n" << obj_vec_.transpose() << std::endl;
  ofs << "eq_mat:\n" << eq_mat_ << std::endl;
  ofs << "eq_vec:\n" << eq_vec_.transpose() << std::endl;
  ofs << "ineq_mat:\n" << ineq_mat_ << std::endl;
  ofs << "ineq_vec:\n" << ineq_vec_.transpose() << std::endl;
  ofs << "x_min:\n" << x_min_.transpose() << std::endl;
  ofs << "x_max:\n" << x_max_.transpose() << std::endl;
}

void QpSolver::printInfo(bool, const std::string & header) const
{
  QSC_INFO_STREAM(header << "QP solver: " << std::to_string(type_));
}

Eigen::VectorXd QpSolver::solve(QpCoeff & qp_coeff)
{
  return solve(qp_coeff.dim_var_, qp_coeff.dim_eq_, qp_coeff.dim_ineq_, qp_coeff.obj_mat_, qp_coeff.obj_vec_,
               qp_coeff.eq_mat_, qp_coeff.eq_vec_, qp_coeff.ineq_mat_, qp_coeff.ineq_vec_, qp_coeff.x_min_,
               qp_coeff.x_max_);
}

#if ENABLE_QLD
QpSolverQld::QpSolverQld()
{
  type_ = QpSolverType::QLD;
  qld_ = std::make_unique<Eigen::QLDDirect>();
}

Eigen::VectorXd QpSolverQld::solve(int dim_var,
                                   int dim_eq,
                                   int dim_ineq,
                                   Eigen::Ref<Eigen::MatrixXd> Q,
                                   const Eigen::Ref<const Eigen::VectorXd> & c,
                                   const Eigen::Ref<const Eigen::MatrixXd> & A,
                                   const Eigen::Ref<const Eigen::VectorXd> & b,
                                   const Eigen::Ref<const Eigen::MatrixXd> & C,
                                   const Eigen::Ref<const Eigen::VectorXd> & d,
                                   const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                   const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  Eigen::MatrixXd AC(dim_eq + dim_ineq, dim_var);
  Eigen::VectorXd bd(dim_eq + dim_ineq);
  AC << -A, -C;
  bd << b, d;

  qld_->problem(dim_var, dim_eq, dim_ineq);
  qld_->solve(Q, c, AC, bd, x_min, x_max, dim_eq);

  if(qld_->fail() == 0)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    QSC_WARN_STREAM("[QpSolverQld::solve] Failed to solve: " << qld_->fail());
  }

  return qld_->result();
}
#endif

#if ENABLE_QUADPROG
QpSolverQuadprog::QpSolverQuadprog()
{
  type_ = QpSolverType::QuadProg;
  quadprog_ = std::make_unique<Eigen::QuadProgDense>();
}

Eigen::VectorXd QpSolverQuadprog::solve(int dim_var,
                                        int dim_eq,
                                        int dim_ineq,
                                        Eigen::Ref<Eigen::MatrixXd> Q,
                                        const Eigen::Ref<const Eigen::VectorXd> & c,
                                        const Eigen::Ref<const Eigen::MatrixXd> & A,
                                        const Eigen::Ref<const Eigen::VectorXd> & b,
                                        const Eigen::Ref<const Eigen::MatrixXd> & C,
                                        const Eigen::Ref<const Eigen::VectorXd> & d,
                                        const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                        const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  int dim_ineq_with_bound = dim_ineq + 2 * dim_var;
  Eigen::MatrixXd C_with_bound(dim_ineq_with_bound, dim_var);
  Eigen::VectorXd d_with_bound(dim_ineq_with_bound);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim_var, dim_var);
  C_with_bound << C, I, -I;
  d_with_bound << d, x_max, -x_min;

  quadprog_->problem(dim_var, dim_eq, dim_ineq_with_bound);
  quadprog_->solve(Q, c, A, b, C_with_bound, d_with_bound);

  if(quadprog_->fail() == 0)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    ROS_WARN_STREAM("[QpSolverQuadprog::solve] Failed to solve: " << quadprog_->fail());
  }

  return quadprog_->result();
}
#endif

#if ENABLE_LSSOL
QpSolverLssol::QpSolverLssol()
{
  type_ = QpSolverType::LSSOL;
  lssol_ = std::make_unique<Eigen::LSSOL_QP>();
}

Eigen::VectorXd QpSolverLssol::solve(int dim_var,
                                     int dim_eq,
                                     int dim_ineq,
                                     Eigen::Ref<Eigen::MatrixXd> Q,
                                     const Eigen::Ref<const Eigen::VectorXd> & c,
                                     const Eigen::Ref<const Eigen::MatrixXd> & A,
                                     const Eigen::Ref<const Eigen::VectorXd> & b,
                                     const Eigen::Ref<const Eigen::MatrixXd> & C,
                                     const Eigen::Ref<const Eigen::VectorXd> & d,
                                     const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                     const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  Eigen::MatrixXd AC(dim_eq + dim_ineq, dim_var);
  Eigen::VectorXd bd_min(dim_eq + dim_ineq);
  Eigen::VectorXd bd_max(dim_eq + dim_ineq);
  AC << A, C;
  bd_min << b, Eigen::VectorXd::Constant(dim_ineq, -1 * std::numeric_limits<double>::infinity());
  bd_max << b, d;

  lssol_->resize(dim_var, dim_eq + dim_ineq, Eigen::lssol::QP2);

  lssol_->persistence(!solve_failed_);
  lssol_->warm(!solve_failed_);

  lssol_->solve(x_min, x_max, Q, c, AC, bd_min, bd_max);

  if(lssol_->inform() == Eigen::lssol::STRONG_MINIMUM)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    std::stringstream sstream;
    lssol_->inform(sstream);
    QSC_WARN_STREAM("[QpSolverLssol::solve] Failed to solve: " << sstream.str());
  }

  return lssol_->result();
}
#endif

#if ENABLE_JRLQP
QpSolverJrlqp::QpSolverJrlqp()
{
  type_ = QpSolverType::JRLQP;
  jrlqp_ = std::make_unique<jrl::qp::GoldfarbIdnaniSolver>();
}

Eigen::VectorXd QpSolverJrlqp::solve(int dim_var,
                                     int dim_eq,
                                     int dim_ineq,
                                     Eigen::Ref<Eigen::MatrixXd> Q,
                                     const Eigen::Ref<const Eigen::VectorXd> & c,
                                     const Eigen::Ref<const Eigen::MatrixXd> & A,
                                     const Eigen::Ref<const Eigen::VectorXd> & b,
                                     const Eigen::Ref<const Eigen::MatrixXd> & C,
                                     const Eigen::Ref<const Eigen::VectorXd> & d,
                                     const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                     const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  Eigen::MatrixXd AC(dim_eq + dim_ineq, dim_var);
  Eigen::VectorXd bd_min(dim_eq + dim_ineq);
  Eigen::VectorXd bd_max(dim_eq + dim_ineq);
  AC << A, C;
  bd_min << b, Eigen::VectorXd::Constant(dim_ineq, -1 * std::numeric_limits<double>::infinity());
  bd_max << b, d;

  jrlqp_->resize(dim_var, dim_eq + dim_ineq, true);

  if(solve_failed_)
  {
    jrlqp_->resetActiveSet();
  }
  else
  {
    jrl::qp::SolverOptions solver_option;
    solver_option.warmStart(true);
    jrlqp_->options(solver_option);
  }

  jrl::qp::TerminationStatus status = jrlqp_->solve(Q, c, AC.transpose(), bd_min, bd_max, x_min, x_max);

  if(status == jrl::qp::TerminationStatus::SUCCESS)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    QSC_WARN_STREAM("[QpSolverJrlqp::solve] Failed to solve: " << status);
  }

  return jrlqp_->solution();
}
#endif

#if ENABLE_QPOASES
QpSolverQpoases::QpSolverQpoases()
{
  type_ = QpSolverType::qpOASES;
}

Eigen::VectorXd QpSolverQpoases::solve(int dim_var,
                                       int dim_eq,
                                       int dim_ineq,
                                       Eigen::Ref<Eigen::MatrixXd> Q,
                                       const Eigen::Ref<const Eigen::VectorXd> & c,
                                       const Eigen::Ref<const Eigen::MatrixXd> & A,
                                       const Eigen::Ref<const Eigen::VectorXd> & b,
                                       const Eigen::Ref<const Eigen::MatrixXd> & C,
                                       const Eigen::Ref<const Eigen::VectorXd> & d,
                                       const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                       const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> AC_row_major(dim_eq + dim_ineq, dim_var);
  Eigen::VectorXd bd_min(dim_eq + dim_ineq);
  Eigen::VectorXd bd_max(dim_eq + dim_ineq);
  AC_row_major << A, C;
  bd_min << b, Eigen::VectorXd::Constant(dim_ineq, -1 * std::numeric_limits<double>::infinity());
  bd_max << b, d;

  qpOASES::returnValue status = qpOASES::TERMINAL_LIST_ELEMENT;
  if(!solve_failed_ && !force_initialize_ && qpoases_ && qpoases_->getNV() == dim_var
     && qpoases_->getNC() == dim_eq + dim_ineq)
  {
    status = qpoases_->hotstart(
        // Since Q is a symmetric matrix, row/column-majors are interchangeable
        Q.data(), c.data(), AC_row_major.data(), x_min.data(), x_max.data(), bd_min.data(), bd_max.data(), n_wsr_);
  }
  if(status != qpOASES::SUCCESSFUL_RETURN)
  {
    qpoases_ = std::make_unique<qpOASES::SQProblem>(dim_var, dim_eq + dim_ineq);
    qpoases_->setPrintLevel(qpOASES::PL_LOW);
    status = qpoases_->init(
        // Since Q is a symmetric matrix, row/column-majors are interchangeable
        Q.data(), c.data(), AC_row_major.data(), x_min.data(), x_max.data(), bd_min.data(), bd_max.data(), n_wsr_);
  }

  if(status == qpOASES::SUCCESSFUL_RETURN)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    QSC_WARN_STREAM("[QpSolverQpoases::solve] Failed to solve: " << static_cast<int>(status));
  }

  Eigen::VectorXd sol(dim_var);
  qpoases_->getPrimalSolution(sol.data());
  return sol;
}
#endif

#if ENABLE_OSQP
QpSolverOsqp::QpSolverOsqp()
{
  type_ = QpSolverType::OSQP;
  osqp_ = std::make_unique<OsqpEigen::Solver>();
}

Eigen::VectorXd QpSolverOsqp::solve(int dim_var,
                                    int dim_eq,
                                    int dim_ineq,
                                    Eigen::Ref<Eigen::MatrixXd> Q,
                                    const Eigen::Ref<const Eigen::VectorXd> & c,
                                    const Eigen::Ref<const Eigen::MatrixXd> & A,
                                    const Eigen::Ref<const Eigen::VectorXd> & b,
                                    const Eigen::Ref<const Eigen::MatrixXd> & C,
                                    const Eigen::Ref<const Eigen::VectorXd> & d,
                                    const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                    const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  int dim_eq_ineq_with_bound = dim_eq + dim_ineq + dim_var;
  Eigen::MatrixXd AC_with_bound(dim_eq_ineq_with_bound, dim_var);
  Eigen::VectorXd bd_with_bound_min(dim_eq_ineq_with_bound);
  Eigen::VectorXd bd_with_bound_max(dim_eq_ineq_with_bound);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim_var, dim_var);
  AC_with_bound << A, C, I;
  bd_with_bound_min << b, Eigen::VectorXd::Constant(dim_ineq, -1 * std::numeric_limits<double>::infinity()), x_min;
  bd_with_bound_max << b, d, x_max;

  auto sparse_start_time = clock::now();
  // Matrices and vectors must be hold during solver's lifetime
  Q_sparse_ = Q.sparseView();
  AC_with_bound_sparse_ = AC_with_bound.sparseView();
  // You must pass unconst vectors to OSQP
  c_ = c;
  bd_with_bound_min_ = bd_with_bound_min;
  bd_with_bound_max_ = bd_with_bound_max;
  auto sparse_end_time = clock::now();
  sparse_duration_ =
      1e3 * std::chrono::duration_cast<std::chrono::duration<double>>(sparse_end_time - sparse_start_time).count();

  // osqp_->settings()->setAbsoluteTolerance(1e-2);
  // osqp_->settings()->setRelativeTolerance(1e-2);
  // osqp_->settings()->setScaledTerimination(1);
  // QSC_INFO_STREAM("max_iter: " << osqp_->settings()->getSettings()->max_iter << ", " <<
  //                 "eps_abs: " << osqp_->settings()->getSettings()->eps_abs << ", " <<
  //                 "eps_rel: " << osqp_->settings()->getSettings()->eps_rel << ", " <<
  //                 "scaled_termination: " << osqp_->settings()->getSettings()->scaled_termination);

  osqp_->settings()->setVerbosity(false);
  osqp_->settings()->setWarmStart(true);
  if(!solve_failed_ && !force_initialize_ && osqp_->isInitialized() && dim_var == osqp_->data()->getData()->n
     && dim_eq_ineq_with_bound == osqp_->data()->getData()->m)
  {
    // Update only matrices and vectors
    osqp_->updateHessianMatrix(Q_sparse_);
    osqp_->updateGradient(c_);
    osqp_->updateLinearConstraintsMatrix(AC_with_bound_sparse_);
    osqp_->updateBounds(bd_with_bound_min_, bd_with_bound_max_);
  }
  else
  {
    // Initialize fully
    if(osqp_->isInitialized())
    {
      osqp_->clearSolver();
      osqp_->data()->clearHessianMatrix();
      osqp_->data()->clearLinearConstraintsMatrix();
    }

    osqp_->data()->setNumberOfVariables(dim_var);
    osqp_->data()->setNumberOfConstraints(dim_eq_ineq_with_bound);
    osqp_->data()->setHessianMatrix(Q_sparse_);
    osqp_->data()->setGradient(c_);
    osqp_->data()->setLinearConstraintsMatrix(AC_with_bound_sparse_);
    osqp_->data()->setLowerBound(bd_with_bound_min_);
    osqp_->data()->setUpperBound(bd_with_bound_max_);
    osqp_->initSolver();
  }

  osqp_->solve();

  if(osqp_->workspace()->info->status_val == OSQP_SOLVED)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    QSC_WARN_STREAM("[QpSolverOsqp::solve] Failed to solve: " << osqp_->workspace()->info->status);
  }

  return osqp_->getSolution();
}
#endif

#if ENABLE_NASOQ
QpSolverNasoq::QpSolverNasoq()
{
  type_ = QpSolverType::NASOQ;
}

Eigen::VectorXd QpSolverNasoq::solve(int dim_var,
                                     int dim_eq,
                                     int dim_ineq,
                                     Eigen::Ref<Eigen::MatrixXd> Q,
                                     const Eigen::Ref<const Eigen::VectorXd> & c,
                                     const Eigen::Ref<const Eigen::MatrixXd> & A,
                                     const Eigen::Ref<const Eigen::VectorXd> & b,
                                     const Eigen::Ref<const Eigen::MatrixXd> & C,
                                     const Eigen::Ref<const Eigen::VectorXd> & d,
                                     const Eigen::Ref<const Eigen::VectorXd> & x_min,
                                     const Eigen::Ref<const Eigen::VectorXd> & x_max)
{
  int dim_ineq_with_bound = dim_ineq + 2 * dim_var;
  Eigen::MatrixXd C_with_bound(dim_ineq_with_bound, dim_var);
  Eigen::VectorXd d_with_bound(dim_ineq_with_bound);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim_var, dim_var);
  C_with_bound << C, I, -I;
  d_with_bound << d, x_max, -x_min;

  auto sparse_start_time = clock::now();
  // Matrices and vectors must be hold during solver's lifetime
  Q_sparse_ = Q.sparseView();
  A_sparse_ = A.sparseView();
  C_with_bound_sparse_ = C_with_bound.sparseView();
  auto sparse_end_time = clock::now();
  sparse_duration_ =
      1e3 * std::chrono::duration_cast<std::chrono::duration<double>>(sparse_end_time - sparse_start_time).count();

  Eigen::VectorXd sol(dim_var), dual_eq(dim_eq), dual_ineq(dim_ineq_with_bound);
  nasoq::QPSettings settings;
  int solve_ret = nasoq::quadprog(Q_sparse_.triangularView<Eigen::Lower>(), c, A_sparse_, b, C_with_bound_sparse_,
                                  d_with_bound, sol, dual_eq, dual_ineq, &settings);

  if(solve_ret == nasoq::Optimal)
  {
    solve_failed_ = false;
  }
  else
  {
    solve_failed_ = true;
    QSC_WARN_STREAM("[QpSolverNasoq::solve] Failed to solve: " << solve_ret);
  }

  return sol;
}
#endif

QpSolverType QpSolverCollection::getAnyQpSolverType()
{
  if(ENABLE_QLD)
  {
    return QpSolverType::QLD;
  }
  else if(ENABLE_QUADPROG)
  {
    return QpSolverType::QuadProg;
  }
  else if(ENABLE_LSSOL)
  {
    return QpSolverType::LSSOL;
  }
  else if(ENABLE_JRLQP)
  {
    return QpSolverType::JRLQP;
  }
  else if(ENABLE_QPOASES)
  {
    return QpSolverType::qpOASES;
  }
  else if(ENABLE_OSQP)
  {
    return QpSolverType::OSQP;
  }
  else if(ENABLE_NASOQ)
  {
    return QpSolverType::NASOQ;
  }
  else
  {
    throw std::runtime_error("[getAnyQpSolverType] No QP solver is enabled.");
    return QpSolverType::Uninitialized;
  }
}

bool QpSolverCollection::isQpSolverEnabled(const QpSolverType & qp_solver_type)
{
  if(qp_solver_type == QpSolverType::Any)
  {
    return getAnyQpSolverType() != QpSolverType::Uninitialized;
  }
  else if(qp_solver_type == QpSolverType::QLD)
  {
    return ENABLE_QLD;
  }
  else if(qp_solver_type == QpSolverType::QuadProg)
  {
    return ENABLE_QUADPROG;
  }
  else if(qp_solver_type == QpSolverType::LSSOL)
  {
    return ENABLE_LSSOL;
  }
  else if(qp_solver_type == QpSolverType::JRLQP)
  {
    return ENABLE_JRLQP;
  }
  else if(qp_solver_type == QpSolverType::qpOASES)
  {
    return ENABLE_QPOASES;
  }
  else if(qp_solver_type == QpSolverType::OSQP)
  {
    return ENABLE_OSQP;
  }
  else if(qp_solver_type == QpSolverType::NASOQ)
  {
    return ENABLE_NASOQ;
  }
  else
  {
    QSC_ERROR_STREAM("[isQpSolverEnabled] Unsupported QP solver: " << std::to_string(static_cast<int>(qp_solver_type)));
    return false;
  }
}

std::shared_ptr<QpSolver> QpSolverCollection::allocateQpSolver(const QpSolverType & qp_solver_type)
{
  std::shared_ptr<QpSolver> qp;

  if(qp_solver_type == QpSolverType::Any)
  {
    return allocateQpSolver(getAnyQpSolverType());
  }
  else if(qp_solver_type == QpSolverType::QLD)
  {
#if ENABLE_QLD
    qp = std::make_shared<QpSolverQld>();
#endif
  }
  else if(qp_solver_type == QpSolverType::QuadProg)
  {
#if ENABLE_QUADPROG
    qp = std::make_shared<QpSolverQuadprog>();
#endif
  }
  else if(qp_solver_type == QpSolverType::LSSOL)
  {
#if ENABLE_LSSOL
    qp = std::make_shared<QpSolverLssol>();
#endif
  }
  else if(qp_solver_type == QpSolverType::JRLQP)
  {
#if ENABLE_JRLQP
    qp = std::make_shared<QpSolverJrlqp>();
#endif
  }
  else if(qp_solver_type == QpSolverType::qpOASES)
  {
#if ENABLE_QPOASES
    qp = std::make_shared<QpSolverQpoases>();
#endif
  }
  else if(qp_solver_type == QpSolverType::OSQP)
  {
#if ENABLE_OSQP
    qp = std::make_shared<QpSolverOsqp>();
#endif
  }
  else if(qp_solver_type == QpSolverType::NASOQ)
  {
#if ENABLE_NASOQ
    qp = std::make_shared<QpSolverNasoq>();
#endif
  }
  else
  {
    QSC_ERROR_STREAM("[allocateQpSolver] Unsupported QP solver: " << std::to_string(static_cast<int>(qp_solver_type)));
  }

  if(!qp)
  {
    QSC_ERROR_STREAM("[allocateQpSolver] Failed to initialize QP solver: " << std::to_string(qp_solver_type));
  }

  return qp;
}
