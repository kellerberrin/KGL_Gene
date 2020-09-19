//
// Created by kellerberrin on 18/9/20.
//

#include "kel_optimize.h"
#include "kel_exec_env.h"

// Must be linked against libnlopt
#include <nlopt.hpp>

namespace kel = kellerberrin;



void kel::Optimize::addEqualityNonLinearConstraint(NonLinearConstraintFunction constraint_function, const std::vector<double>& data, double tolerance) {

  NonLinearConstraint constraint;
  constraint.constraint_function = constraint_function;
  constraint.data = data;
  constraint.tolerance = tolerance;
  equality_constraints_.push_back(constraint);

}


void kel::Optimize::addInequalityNonLinearConstraint(NonLinearConstraintFunction constraint_function, const std::vector<double>& data, double tolerance) {

  NonLinearConstraint constraint;
  constraint.constraint_function = constraint_function;
  constraint.data = data;
  constraint.tolerance = tolerance;
  inequality_constraints_.push_back(constraint);

}


std::pair<kel::OptimizationResult, double> kel::Optimize::optimize(const std::vector<double>& inital_point,
                                                                   std::vector<double>& final_point) {

  auto opt_alg = static_cast<nlopt::algorithm>(convertAlgorithm(opt_alg_));
  nlopt::opt opt(opt_alg, dimension_);

  if (not upper_bound_.empty()) {

    opt.set_upper_bounds(upper_bound_);

  }

  if (not lower_bound_.empty()) {

    opt.set_lower_bounds(lower_bound_);

  }

  // Setup the objective function
  nlopt::vfunc objective = *objective_.target<nlopt::vfunc>();
  if (opt_type_ == OptimizationType::MINIMIZE) {

    opt.set_min_objective(objective, nullptr);

  } else {

    opt.set_max_objective(objective, nullptr);

  }

  // Setup any equality constraints
  for (auto constraint: equality_constraints_) {

    nlopt::vfunc constraint_func = *constraint.constraint_function.target<nlopt::vfunc>();
    opt.add_equality_constraint(constraint_func, &constraint.data.front(), constraint.tolerance);

  }

  // Setup any inequality constraints
  for (auto constraint: inequality_constraints_) {

    nlopt::vfunc constraint_func = *constraint.constraint_function.target<nlopt::vfunc>();
    opt.add_inequality_constraint(constraint_func, &constraint.data.front(), constraint.tolerance);

  }

  opt.set_xtol_rel(1e-4);
  final_point = inital_point;

  try{

    double minf;
    nlopt::result result = opt.optimize(final_point, minf);
    return {static_cast<OptimizationResult>(result), minf};

  }
  catch(std::exception &e) {

    ExecEnv::log().error("Optimize:: nlopt failed: {}", e.what());
    return {OptimizationResult::FAILURE, 0.0};

  }

}

void kel::Optimize::boundingHypercube(const std::vector<double>& upper_bound, const std::vector<double>& lower_bound) {

  upper_bound_ = upper_bound;
  lower_bound_ = lower_bound;

}

size_t kel::Optimize::convertAlgorithm(OptimizationAlgorithm algorithm) {

  switch(algorithm) {

    case OptimizationAlgorithm::GN_DIRECT: return nlopt::GN_DIRECT;
    case OptimizationAlgorithm::GN_DIRECT_L: return nlopt::GN_DIRECT_L;
    case OptimizationAlgorithm::GN_DIRECT_L_RAND: return nlopt::GN_DIRECT_L_RAND;
    case OptimizationAlgorithm::GN_DIRECT_NOSCAL: return nlopt::GN_DIRECT_NOSCAL;
    case OptimizationAlgorithm::GN_DIRECT_L_NOSCAL: return nlopt::GN_DIRECT_L_NOSCAL;
    case OptimizationAlgorithm::GN_DIRECT_L_RAND_NOSCAL: return nlopt::GN_DIRECT_L_RAND_NOSCAL;
    case OptimizationAlgorithm::GN_ORIG_DIRECT: return nlopt::GN_ORIG_DIRECT;
    case OptimizationAlgorithm::GN_ORIG_DIRECT_L: return nlopt::GN_ORIG_DIRECT_L;
    case OptimizationAlgorithm::GD_STOGO: return nlopt::GD_STOGO;
    case OptimizationAlgorithm::GD_STOGO_RAND: return nlopt::GD_STOGO_RAND;
    case OptimizationAlgorithm::LD_LBFGS_NOCEDAL: return nlopt::LD_LBFGS_NOCEDAL;
    case OptimizationAlgorithm::LD_LBFGS: return nlopt::LD_LBFGS;
    case OptimizationAlgorithm::LN_PRAXIS: return nlopt::LN_PRAXIS;
    case OptimizationAlgorithm::LD_VAR1: return nlopt::LD_VAR1;
    case OptimizationAlgorithm::LD_VAR2: return nlopt::LD_VAR2;
    case OptimizationAlgorithm::LD_TNEWTON: return nlopt::LD_TNEWTON;
    case OptimizationAlgorithm::LD_TNEWTON_RESTART: return nlopt::LD_TNEWTON_RESTART;
    case OptimizationAlgorithm::LD_TNEWTON_PRECOND: return nlopt::LD_TNEWTON_PRECOND;
    case OptimizationAlgorithm::LD_TNEWTON_PRECOND_RESTART: return nlopt::LD_TNEWTON_PRECOND_RESTART;
    case OptimizationAlgorithm::GN_CRS2_LM: return nlopt::GN_CRS2_LM;
    case OptimizationAlgorithm::GN_MLSL: return nlopt::GN_MLSL;
    case OptimizationAlgorithm::GD_MLSL: return nlopt::GD_MLSL;
    case OptimizationAlgorithm::GN_MLSL_LDS: return nlopt::GN_MLSL_LDS;
    case OptimizationAlgorithm::GD_MLSL_LDS: return nlopt::GD_MLSL_LDS;
    case OptimizationAlgorithm::LD_MMA: return nlopt::LD_MMA;
    case OptimizationAlgorithm::LN_COBYLA: return nlopt::LN_COBYLA;
    case OptimizationAlgorithm::LN_NEWUOA: return nlopt::LN_NEWUOA;
    case OptimizationAlgorithm::LN_NEWUOA_BOUND: return nlopt::LN_NEWUOA_BOUND;
    case OptimizationAlgorithm::LN_NELDERMEAD: return nlopt::LN_NELDERMEAD;
    case OptimizationAlgorithm::LN_SBPLX: return nlopt::LN_SBPLX;
    case OptimizationAlgorithm::LN_AUGLAG: return nlopt::LN_AUGLAG;
    case OptimizationAlgorithm::LD_AUGLAG: return nlopt::LD_AUGLAG;
    case OptimizationAlgorithm::LN_AUGLAG_EQ: return nlopt::LN_AUGLAG_EQ;
    case OptimizationAlgorithm::LD_AUGLAG_EQ: return nlopt::LD_AUGLAG_EQ;
    case OptimizationAlgorithm::LN_BOBYQA: return nlopt::LN_BOBYQA;
    case OptimizationAlgorithm::GN_ISRES: return nlopt::GN_ISRES;
    case OptimizationAlgorithm::AUGLAG: return nlopt::AUGLAG;
    case OptimizationAlgorithm::AUGLAG_EQ: return nlopt::AUGLAG_EQ;
    case OptimizationAlgorithm::G_MLSL: return nlopt::G_MLSL;
    case OptimizationAlgorithm::G_MLSL_LDS: return nlopt::G_MLSL_LDS;
    case OptimizationAlgorithm::LD_SLSQP: return nlopt::LD_SLSQP;
    case OptimizationAlgorithm::LD_CCSAQ: return nlopt::LD_CCSAQ;
    case OptimizationAlgorithm::GN_ESCH: return nlopt::GN_ESCH;
    case OptimizationAlgorithm::GN_AGS: return nlopt::GN_AGS;

  }

  ExecEnv::log().error("Optimize::convertAlgorithm; Invalid Optimization Algorithm specified");
  return nlopt::NUM_ALGORITHMS;

}

std::string kel::Optimize::returnDescription(OptimizationResult result) {

  switch (result) {

    case OptimizationResult::FAILURE: return "FAILURE";       /* generic failure code */
    case OptimizationResult::INVALID_ARGS: return "INVALID_ARGS";
    case OptimizationResult::OUT_OF_MEMORY: return "OUT_OF_MEMORY";
    case OptimizationResult::ROUNDOFF_LIMITED: return "ROUNDOFF_LIMITED";
    case OptimizationResult::FORCED_STOP: return "FORCED_STOP";
    case OptimizationResult::SUCCESS: return "SUCCESS";           /* generic success code */
    case OptimizationResult::STOPVAL_REACHED: return "STOPVAL_REACHED";
    case OptimizationResult::FTOL_REACHED: return "FTOL_REACHED";
    case OptimizationResult::XTOL_REACHED: return "XTOL_REACHED";
    case OptimizationResult::MAXEVAL_REACHED: return "MAXEVAL_REACHED";
    case OptimizationResult::MAXTIME_REACHED: return "MAXTIME_REACHED";

  }

  return "UNKNOWN - INVALID RETURNCODE";

}
