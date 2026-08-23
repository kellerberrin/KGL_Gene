// Kellerberrin 2026.

#include "kel_optimize.h"

#include <cmath>
#include <limits>

// Must be linked against "libnlopt"
#include <nlopt.hpp>

namespace kel = kellerberrin;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local test functions (NLopt example problem).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {

double myvconstraint(const std::vector<double>& x, std::vector<double>& grad, std::vector<double>& data) {
  const double a = data[0], b = data[1];
  if (!grad.empty()) {
    grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
    grad[1] = -1.0;
  }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

double myvfunc(const std::vector<double>& x, std::vector<double>& grad) {
  if (!grad.empty()) {
    grad[0] = 0.0;
    grad[1] = 0.5 / std::sqrt(x[1]);
  }
  return std::sqrt(x[1]);
}

} // anonymous namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Experimental test optimization (to be removed).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kel::Optimize::opt_test() {

  std::vector<double> x_values{1.234, 5.678};
  Optimize opt(OptimizationAlgorithm::LD_MMA, x_values.size(), OptimizationType::MINIMIZE);
  std::vector<double> lower_bound{std::numeric_limits<double>::lowest(), 0};
  opt.boundingHypercube({}, lower_bound);
  opt.addInequalityNonLinearConstraint(myvconstraint, std::vector<double>{2, 0}, 1e-8);
  opt.addInequalityNonLinearConstraint(myvconstraint, std::vector<double>{-1, 1}, 1e-8);
  opt.stoppingCriteria(OptimizeStoppingType::RELATIVE_PARAMETER_THRESHOLD, {1e-04});
  OptDerivObjectiveFn obj = myvfunc;
  auto [result_code, optimal_value, iterations] = opt.optimize(x_values, obj);
  ExecEnv::log().info("Optimize::*****************************************************************************************************");
  ExecEnv::log().info("Optimize:: Found minimum at f({},{}) = {}, optimizer result: {}, iterations: {}",
      x_values[0], x_values[1], optimal_value, Optimize::returnDescription(result_code), iterations);
  std::vector<double> x_values2{5.678, 1.234};
  auto [result_code2, optimal_value2, iterations2] = opt.optimize(x_values2, myvfunc);
  ExecEnv::log().info("Optimize:: Found minimum at f({},{}) = {}, optimizer result: {}, iterations: {}",
      x_values2[0], x_values2[1], optimal_value2, Optimize::returnDescription(result_code2), iterations2);
  ExecEnv::log().info("Optimize::*****************************************************************************************************");

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Equality constraints (no-data variants)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kel::Optimize::addEqualityNonLinearConstraint(OptDerivConstraintFn constraint_function, double tolerance) {

  ConstraintFunction closure = [constraint_function = std::move(constraint_function)](const std::vector<double>& x, std::vector<double>& grad) -> double {
    return constraint_function(x, grad);
  };
  equality_constraints_.push_back({std::move(closure), tolerance});

}

void kel::Optimize::addEqualityNonLinearConstraint(OptConstraintFn constraint_function, double tolerance) {

  ConstraintFunction closure = [constraint_function = std::move(constraint_function)](const std::vector<double>& x, std::vector<double>& grad) -> double {
    if (not grad.empty()) {

      ExecEnv::log().error("Optimize::addEqualityNonLinearConstraint; expected non-differential algorithm");
      return std::numeric_limits<double>::quiet_NaN();

    }
    return constraint_function(x);
  };
  equality_constraints_.push_back({std::move(closure), tolerance});

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Inequality constraints (no-data variants)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kel::Optimize::addInequalityNonLinearConstraint(OptDerivConstraintFn constraint_function, double tolerance) {

  ConstraintFunction closure = [constraint_function = std::move(constraint_function)](const std::vector<double>& x, std::vector<double>& grad) -> double {
    return constraint_function(x, grad);
  };
  inequality_constraints_.push_back({std::move(closure), tolerance});

}

void kel::Optimize::addInequalityNonLinearConstraint(OptConstraintFn constraint_function, double tolerance) {

  ConstraintFunction closure = [constraint_function = std::move(constraint_function)](const std::vector<double>& x, std::vector<double>& grad) -> double {
    if (not grad.empty()) {

      ExecEnv::log().error("Optimize::addInequalityNonLinearConstraint; expected non-differential algorithm");
      return std::numeric_limits<double>::quiet_NaN();

    }
    return constraint_function(x);
  };
  inequality_constraints_.push_back({std::move(closure), tolerance});

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stopping criteria — now uses AND (`&&`) for dimension checks
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kel::Optimize::stoppingCriteria(OptimizeStoppingType stopping_type, const std::vector<double>& stopping_value) {

  const auto require_size = [&](std::size_t expected) {
    if (stopping_value.size() != expected) {

      ExecEnv::log().error("Optimize::stoppingCriteria; stopping criteria requires vector size {}, actual size: {}",
                           expected, stopping_value.size());
      return false;

    }
    return true;
  };

  const auto require_size1_or_dim = [&]() {
    if (stopping_value.size() != 1 && stopping_value.size() != dimension_) {

      ExecEnv::log().error("Optimize::stoppingCriteria; stopping criteria requires vector size 1 OR objective function dimension: {}, actual size: {}",
                           dimension_, stopping_value.size());
      return false;

    }
    return true;
  };

  switch(stopping_type) {

    case OptimizeStoppingType::FUNCTION_VALUE:               // Stop when a particular function value is reached.
    case OptimizeStoppingType::RELATIVE_FUNCTION_THRESHOLD:   // Stop when the relative objective function update is below a threshold.
    case OptimizeStoppingType::ABSOLUTE_FUNCTION_THRESHOLD:    // Stop when the absolute objective function update is below a threshold.
    case OptimizeStoppingType::RELATIVE_PARAMETER_THRESHOLD:  // Stop when the relative weighted (normed) parameter vector update is below a threshold.
    case OptimizeStoppingType::MAXIMUM_EVALUATIONS:           // Stop when the maximum number of evaluations have been reached.
    case OptimizeStoppingType::MAXIMUM_TIME:                  // Stop when the specified time in seconds has elapsed.
    {
      if (!require_size(1)) { return; }
      break;
    }

    case OptimizeStoppingType::RELATIVE_PARAMETER_WEIGHTS:  // Set the relative parameter norm weights.
    case OptimizeStoppingType::ABSOLUTE_PARAMETER_THRESHOLD:  // Stop when the parameter vector update is below a threshold.
    {
      if (!require_size1_or_dim()) { return; }
      break;
    }

  }

  stopping_vector_.push_back({stopping_type, stopping_value});

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Optimize entry points
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

kel::OptResultTuple kel::Optimize::optimize(std::vector<double>& x_parameter_vector, OptDerivObjectiveFn objective) {

  ObjectiveFunction closure = [objective = std::move(objective)](const std::vector<double>& x, std::vector<double>& grad) -> double {
    if (grad.empty()) {

      ExecEnv::log().error("Optimize::optimize; expected differential optimization algorithm");
      return std::numeric_limits<double>::quiet_NaN();

    }
    return objective(x, grad);
  };
  return runOptimize(x_parameter_vector, std::move(closure));

}


kel::OptResultTuple kel::Optimize::optimize(std::vector<double>& x_parameter_vector, OptObjectiveFn objective) {

  ObjectiveFunction closure = [objective = std::move(objective)](const std::vector<double>& x, std::vector<double>& grad) -> double {
    if (not grad.empty()) {

      ExecEnv::log().error("Optimize::optimize; expected non-differential optimization algorithm");
      return std::numeric_limits<double>::quiet_NaN();

    }
    return objective(x);
  };
  return runOptimize(x_parameter_vector, std::move(closure));

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Private entry to the underlying optimizer code.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

kel::OptResultTuple kel::Optimize::runOptimize(std::vector<double>& parameter_x_vector, ObjectiveFunction objective) {

  if (parameter_x_vector.size() != dimension_) {

    ExecEnv::log().error("Optimize::optimize; Objective Function Dimension: {}, Initial point Dimension: {}",
                         dimension_, parameter_x_vector.size());
    return {OptimizationResult::FAILURE, 0.0, 0};

  }

  objective_ = std::move(objective);

  try {

    const auto opt_alg = static_cast<nlopt::algorithm>(convertAlgorithm(opt_alg_));
    nlopt::opt opt(opt_alg, static_cast<unsigned>(dimension_));

    if (not upper_bound_.empty()) {

      opt.set_upper_bounds(upper_bound_);

    }

    if (not lower_bound_.empty()) {

      opt.set_lower_bounds(lower_bound_);

    }

    // Setup the objective function. The trampoline receives `this` as void*.
    auto nlopt_objective = &Optimize::objectiveCallback;
    void* self_ptr = this;
    if (opt_type_ == OptimizationType::MINIMIZE) {

      opt.set_min_objective(nlopt_objective, self_ptr);

    } else {

      opt.set_max_objective(nlopt_objective, self_ptr);

    }

    nlopt::vfunc constraint_func = &Optimize::constraintCallback;
    // Setup any equality constraints
    for (auto& constraint : equality_constraints_) {

      opt.add_equality_constraint(constraint_func, &constraint, constraint.tolerance);

    }

    // Setup any inequality constraints
    for (auto& constraint : inequality_constraints_) {

      opt.add_inequality_constraint(constraint_func, &constraint, constraint.tolerance);

    }

    for (auto& [stopping_type, stopping_value] : stopping_vector_) {

      if (stopping_value.size() != 1 && stopping_value.size() != dimension_) {

        ExecEnv::log().error("Optimize::optimize; Invalid stopping value size: {}", stopping_value.size());
        continue;

      }

      switch(stopping_type) {

        case OptimizeStoppingType::FUNCTION_VALUE:     // Stop when a particular function value is reached.
          opt.set_stopval(stopping_value.front());
          break;

        case OptimizeStoppingType::RELATIVE_FUNCTION_THRESHOLD:   // Stop when the relative objective function update is below a threshold.
          opt.set_ftol_rel(stopping_value.front());
          break;

        case OptimizeStoppingType::ABSOLUTE_FUNCTION_THRESHOLD:    // Stop when the absolute objective function update is below a threshold.
          opt.set_ftol_abs(stopping_value.front());
          break;

        case OptimizeStoppingType::RELATIVE_PARAMETER_THRESHOLD:  // Stop when the relative weighted (normed) parameter vector update is below a threshold.
          opt.set_xtol_rel(stopping_value.front());
          break;

        case OptimizeStoppingType::RELATIVE_PARAMETER_WEIGHTS:  // Set the relative parameter weights.
          if (stopping_value.size() == 1) {

            opt.set_x_weights(stopping_value.front());

          } else {

            opt.set_x_weights(stopping_value);

          }
          break;

        case OptimizeStoppingType::ABSOLUTE_PARAMETER_THRESHOLD:  // Stop when the absolute weighted (normed) parameter vector update is below a threshold.
          if (stopping_value.size() == 1) {

            opt.set_xtol_abs(stopping_value.front());

          } else {

            opt.set_xtol_abs(stopping_value);

          }
          break;

        case OptimizeStoppingType::MAXIMUM_EVALUATIONS: // Stop when the maximum number of evaluations have been reached.
          opt.set_maxeval(static_cast<int>(stopping_value.front()));
          break;

        case OptimizeStoppingType::MAXIMUM_TIME:    // Stop when the specified time in seconds has elapsed.
          opt.set_maxtime(stopping_value.front());
          break;

      }

    }

    std::size_t iterations = 0;

    double optimal_function_value;
    nlopt::result result = opt.optimize(parameter_x_vector, optimal_function_value);
    iterations = static_cast<std::size_t>(opt.get_numevals());
    return {static_cast<OptimizationResult>(result), optimal_function_value, iterations};

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("Optimize:: nlopt failed: {}", e.what());
    return {OptimizationResult::FAILURE, 0.0, 0};

  }

}

void kel::Optimize::boundingHypercube(const std::vector<double>& upper_bound, const std::vector<double>& lower_bound) {

  upper_bound_ = upper_bound;
  lower_bound_ = lower_bound;

}

std::size_t kel::Optimize::convertAlgorithm(OptimizationAlgorithm algorithm) {

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