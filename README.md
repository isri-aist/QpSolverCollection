# qp_solver_collection
Unified interface for quadratic programming solvers

## Supported QP solvers
- [QLD](https://github.com/jrl-umi3218/eigen-qld)
- [QuadProg](https://github.com/jrl-umi3218/eigen-quadprog)
- [LSSOL](https://gite.lirmm.fr/multi-contact/eigen-lssol) (private)
- [JRLQP](https://github.com/jrl-umi3218/jrl-qp)
- [qpOASES](https://github.com/coin-or/qpOASES)
- [OSQP](https://osqp.org/)
- [NASOQ](https://nasoq.github.io/)

## Installation

### Installation procedure

None of the following QP solvers are mandatory to install. You can select the desired QP solvers.

#### QLD
Install [eigen-qld](https://github.com/jrl-umi3218/eigen-qld).
Add `-DENABLE_QLD=ON` to the catkin build command.

#### QuadProg
Install [eigen-quadprog](https://github.com/jrl-umi3218/eigen-quadprog).
Add `-DENABLE_QUADPROG=ON` to the catkin build command.

#### LSSOL
Install [eigen-lssol](https://gite.lirmm.fr/multi-contact/eigen-lssol).
Add `-DENABLE_LSSOL=ON` to the catkin build command.

#### JRLQP
Install [BlockStructure branch](https://github.com/jrl-umi3218/jrl-qp/tree/topic/BlockStructure) of jrl-qp.
Add `-DENABLE_JRLQP=ON` to the catkin build command.

#### qpOASES
Install master branch of [qpOASES](https://github.com/coin-or/qpOASES) with `-DBUILD_SHARED_LIBS=ON` cmake option.
Add `-DENABLE_QPOASES=ON` to the catkin build command.
Also, add `-DQPOASES_INCLUDE_DIR=<path to qpOASES.hpp>` and `-DQPOASES_LIBRARY_DIR=<path to libqpOASES.so>` to the catkin build command.

#### OSQP
Install master branch of [osqp](https://github.com/osqp/osqp) and [osqp-eigen](https://github.com/robotology/osqp-eigen).
Add `-DENABLE_OSQP=ON` to the catkin build command.

#### NASOQ
Install [cmake-install branch](https://github.com/mmurooka/nasoq/tree/cmake-install) of nasoq.
Add `-DENABLE_NASOQ=ON` to the catkin build command.
