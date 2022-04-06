# qp_solver_collection
Unified C++ interface for quadratic programming solvers

[![CI](https://github.com/isri-aist/QpSolverCollection/actions/workflows/ci.yaml/badge.svg)](https://github.com/isri-aist/QpSolverCollection/actions/workflows/ci.yaml)
[![Documentation](https://img.shields.io/badge/doxygen-online-brightgreen?logo=read-the-docs&style=flat)](https://isri-aist.github.io/QpSolverCollection/)

## Features
- Unified C++ interface to many QP solvers
- Can be built as a ROS package
- High portability decoupled from each QP solver by [Pimpl technique](https://en.cppreference.com/w/cpp/language/pimpl)

## Supported QP solvers
- [QLD](https://github.com/jrl-umi3218/eigen-qld)
- [QuadProg](https://github.com/jrl-umi3218/eigen-quadprog)
- [JRLQP](https://github.com/jrl-umi3218/jrl-qp)
- [qpOASES](https://github.com/coin-or/qpOASES)
- [OSQP](https://osqp.org/)
- [NASOQ](https://nasoq.github.io/)
- [LSSOL](https://gite.lirmm.fr/multi-contact/eigen-lssol) (private)

## Installation

### Installation procedure
It is assumed that ROS is installed.

1. Install the QP solver you wish to use according to [this section](https://github.com/isri-aist/QpSolverCollection#qp-solver-installation). You can skip installing QP solvers that you do not use.

2. Setup catkin workspace.
```bash
$ mkdir -p ~/ros/ws_qp_solver_collection/src
$ cd ~/ros/ws_qp_solver_collection
$ wstool init src
$ wstool set -t src isri-aist/QpSolverCollection git@github.com:isri-aist/QpSolverCollection.git --git -y
$ wstool update -t src
```

3. Install dependent packages.
```bash
$ source /opt/ros/${ROS_DISTRO}/setup.bash
$ rosdep install -y -r --from-paths src --ignore-src
```

4. Build a package.
```bash
$ catkin build -DCMAKE_BUILD_TYPE=RelWithDebInfo <qp-solver-flags> --catkin-make-args all tests
```
See [this section](https://github.com/isri-aist/QpSolverCollection#qp-solver-installation) for `<qp-solver-flags>`.

### QP solver installation
As all supported QP solvers are installed in [CI](https://github.com/isri-aist/QpSolverCollection/blob/master/.github/workflows/ci.yaml), please refer to the installation procedure.  
Please refer to the license specified in each QP solver when using it.

#### QLD
Install [eigen-qld](https://github.com/jrl-umi3218/eigen-qld).  
Add `-DENABLE_QLD=ON` to the catkin build command (i.e., `<qp-solver-flags>`).

#### QuadProg
Install [eigen-quadprog](https://github.com/jrl-umi3218/eigen-quadprog).  
Add `-DENABLE_QUADPROG=ON` to the catkin build command (i.e., `<qp-solver-flags>`).

#### JRLQP
Install [BlockStructure branch](https://github.com/jrl-umi3218/jrl-qp/tree/topic/BlockStructure) of jrl-qp  .
Add `-DENABLE_JRLQP=ON` to the catkin build command (i.e., `<qp-solver-flags>`).

#### qpOASES
Install master branch of [qpOASES](https://github.com/coin-or/qpOASES) with `-DBUILD_SHARED_LIBS=ON` cmake option.  
Add `-DENABLE_QPOASES=ON` to the catkin build command (i.e., `<qp-solver-flags>`).
Also, add `-DQPOASES_INCLUDE_DIR=<path to qpOASES.hpp>` and `-DQPOASES_LIBRARY_DIR=<path to libqpOASES.so>` to the catkin build command.

#### OSQP
Install master branch of [osqp](https://github.com/osqp/osqp) and [osqp-eigen](https://github.com/robotology/osqp-eigen).  
Add `-DENABLE_OSQP=ON` to the catkin build command (i.e., `<qp-solver-flags>`).

#### NASOQ
Install [cmake-install branch](https://github.com/mmurooka/nasoq/tree/cmake-install) of nasoq.  
Add `-DENABLE_NASOQ=ON` to the catkin build command (i.e., `<qp-solver-flags>`).

#### LSSOL (private)
Install [eigen-lssol](https://gite.lirmm.fr/multi-contact/eigen-lssol).  
Add `-DENABLE_LSSOL=ON` to the catkin build command (i.e., `<qp-solver-flags>`).

## How to use
See [TestSampleQP.cpp](https://github.com/isri-aist/QpSolverCollection/blob/master/tests/TestSampleQP.cpp) for examples of solving QP problems.
