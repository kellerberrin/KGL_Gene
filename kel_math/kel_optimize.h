// Kellerberrin 2026.

#ifndef KEL_OPTIMIZE_H
#define KEL_OPTIMIZE_H

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "kel_exec_env.h"

namespace kellerberrin {   //  organization level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The optimize object is a facade in front of the nlopt optimization library
// If this object is used then the executable must link to "libnlopt".
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Optimization Algorithms. A mix of constrained local and global algorithms, some require
// derivatives to speed convergence. See the nlopt documentation for further description.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class OptimizationAlgorithm {
  GN_DIRECT = 0,
  GN_DIRECT_L,
  GN_DIRECT_L_RAND,
  GN_DIRECT_NOSCAL,
  GN_DIRECT_L_NOSCAL,
  GN_DIRECT_L_RAND_NOSCAL,
  GN_ORIG_DIRECT,
  GN_ORIG_DIRECT_L,
  GD_STOGO,
  GD_STOGO_RAND,
  LD_LBFGS_NOCEDAL,
  LD_LBFGS,
  LN_PRAXIS,
  LD_VAR1,
  LD_VAR2,
  LD_TNEWTON,
  LD_TNEWTON_RESTART,
  LD_TNEWTON_PRECOND,
  LD_TNEWTON_PRECOND_RESTART,
  GN_CRS2_LM,
  GN_MLSL,
  GD_MLSL,
  GN_MLSL_LDS,
  GD_MLSL_LDS,
  LD_MMA,                             // Method of Moving Asymptotes (LD - local, derivative)
  LN_COBYLA,
  LN_NEWUOA,
  LN_NEWUOA_BOUND,
  LN_NELDERMEAD,                    // Nelder-Mead simplex algorithm (LN - local, no-derivative)
  LN_SBPLX,
  LN_AUGLAG,
  LD_AUGLAG,
  LN_AUGLAG_EQ,
  LD_AUGLAG_EQ,
  LN_BOBYQA,
  GN_ISRES,
  AUGLAG,
  AUGLAG_EQ,
  G_MLSL,
  G_MLSL_LDS,
  LD_SLSQP,
  LD_CCSAQ,
  GN_ESCH,
  GN_AGS,
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Optimization return results. Same values as returned by the nlopt package, see documentation.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class OptimizationResult {
  FAILURE = -1,         // generic failure code
  INVALID_ARGS = -2,
  OUT_OF_MEMORY = -3,
  ROUNDOFF_LIMITED = -4,
  FORCED_STOP = -5,
  SUCCESS = 1,          // generic success code
  STOPVAL_REACHED = 2,
  FTOL_REACHED = 3,
  XTOL_REACHED = 4,
  MAXEVAL_REACHED = 5,
  MAXTIME_REACHED = 6,
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// How the objective function is optimized
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class OptimizationType {
  MAXIMIZE = 1,
  MINIMIZE = 2
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Stopping criteria (can be more than 1 criteria)
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class OptimizeStoppingType {
  FUNCTION_VALUE = 1,              // Stop when a particular function value is reached.
  RELATIVE_FUNCTION_THRESHOLD = 2, // Stop when the relative objective function update is below a threshold.
  ABSOLUTE_FUNCTION_THRESHOLD = 3,  // Stop when the absolute objective function update is below a threshold.
  RELATIVE_PARAMETER_THRESHOLD = 4,// Stop when the relative weighted (normed) parameter vector update is below a threshold.
  RELATIVE_PARAMETER_WEIGHTS = 5,  // Set the relative parameter norm weights.
  ABSOLUTE_PARAMETER_THRESHOLD = 6,// Stop when the parameter vector update is below a threshold.
  MAXIMUM_EVALUATIONS = 7,         // Stop when the maximum number of evaluations have been reached.
  MAXIMUM_TIME = 8                 // Stop when the specified time in seconds has elapsed.
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The optimize class, this a thin facade over the corresponding nlopt functionality.
//
// Public API: All user callbacks are std::function objects. Data-carrying variants capture the data
// by value (via shared_ptr) inside a closure, eliminating the void* type-erasure of the original API.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The derivative and data optimization and non-linear constraint function.
template <class OptData> using OptDerivDataObjectiveFn = std::function<double(const std::vector<double>& x, std::vector<double>& grad, OptData& data)>;
template <class OptData> using OptDerivDataConstraintFn = OptDerivDataObjectiveFn<OptData>;

// The no derivative, with data optimization and non-linear constraint function.
template <class OptData> using OptDataObjectiveFn = std::function<double(const std::vector<double>& x, OptData& data)>;
template <class OptData> using OptDataConstraintFn = OptDataObjectiveFn<OptData>;

// The derivative, no data optimization and non-linear constraint function.
using OptDerivObjectiveFn = std::function<double(const std::vector<double>& x, std::vector<double>& grad)>;
using OptDerivConstraintFn = OptDerivObjectiveFn;

// The no derivative, no data optimization and non-linear constraint function.
using OptObjectiveFn = std::function<double(const std::vector<double>& x)>;
using OptConstraintFn = OptObjectiveFn;

// Returned optimization results, the returned tuple is [result_code, optimal_function_value, iterations]
using OptResultTuple = std::tuple<OptimizationResult, double, std::size_t>;

class Optimize {

private:

  // Internal type-erased objective: data is captured inside the std::function closure.
  using ObjectiveFunction = std::function<double(const std::vector<double>& x, std::vector<double>& grad)>;
  using ConstraintFunction = std::function<double(const std::vector<double>& x, std::vector<double>& grad)>;

  struct NonLinearConstraint {
    ConstraintFunction constraint_function;
    double tolerance;
  };

  struct OptimalStopping {
    OptimizeStoppingType stopping_type;
    std::vector<double> stopping_value;
  };

public:

  Optimize(OptimizationAlgorithm opt_alg, std::size_t dimension, OptimizationType opt_type)
    : opt_alg_(opt_alg), dimension_(dimension), opt_type_(opt_type) {}
  ~Optimize() = default;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform the optimizations. The initial parameter vector is updated to the function x parameter stopping value.
  // The returned tuple is [result_code, optimal_function_value, iterations]
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Objective with derivative info, no data.
  OptResultTuple optimize(std::vector<double>& x_parameter_vector, OptDerivObjectiveFn objective);
  // Objective with no derivative, no data.
  OptResultTuple optimize(std::vector<double>& x_parameter_vector, OptObjectiveFn objective);

  // Template function for an optimization function with derivative info and data.
  // The data reference is captured by reference inside the closure; the caller must keep `data` alive
  // for the duration of optimize().
  template<class OptData>
  OptResultTuple optimize(std::vector<double>& x_parameter_vector,
                          OptData& data,
                          OptDerivDataObjectiveFn<OptData> objective) {

    ObjectiveFunction closure = [&data, objective = std::move(objective)](const std::vector<double>& x, std::vector<double>& grad) -> double {
      return objective(x, grad, data);
    };
    return runOptimize(x_parameter_vector, std::move(closure));

  }

  // Template function for an optimization function with data.
  template<class OptData>
  OptResultTuple optimize(std::vector<double>& x_parameter_vector,
                          OptData& data,
                          OptDataObjectiveFn<OptData> objective) {

    ObjectiveFunction closure = [&data, objective = std::move(objective)](const std::vector<double>& x, std::vector<double>& grad) -> double {
      if (not grad.empty()) {

        ExecEnv::log().error("Optimize::optimize<OptData>; expected non-differential optimization algorithm");
        return std::numeric_limits<double>::quiet_NaN();

      }
      return objective(x, data);
    };
    return runOptimize(x_parameter_vector, std::move(closure));

  }

  // Convert the return value to a string
  [[nodiscard]] static std::string returnDescription(OptimizationResult result);
  // Convert the return value to a boolean.
  [[nodiscard]] static bool returnSuccess(OptimizationResult result) noexcept { return static_cast<int>(result) > 0; }

  // Define the hypercube which contains the solution. Must be the same dimension as the objective function.
  // Empty vector is ignored
  void boundingHypercube(const std::vector<double>& upper_bound = {}, const std::vector<double>& lower_bound = {});

  // Stopping Criteria. The stopping criteria vector either has 1 element.
  // Or for x parameter stopping criteria (only), the same dimensionality of the objective function parameters.
  void stoppingCriteria(OptimizeStoppingType stopping_type, const std::vector<double>& stopping_value);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // The non-linear equality functions.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // With derivative info and data. The data is *copied* into the closure so temporaries are allowed.
  template<class OptData>
  void addEqualityNonLinearConstraint(OptDerivDataConstraintFn<OptData> constraint_function,
                                      OptData data,
                                      double tolerance) {

    auto stored = std::make_shared<OptData>(std::move(data));
    ConstraintFunction closure = [constraint_function = std::move(constraint_function), stored](const std::vector<double>& x, std::vector<double>& grad) -> double {
      return constraint_function(x, grad, *stored);
    };
    equality_constraints_.push_back({std::move(closure), tolerance});

  }

  // Without derivative, with data. The data is *copied* into the closure so temporaries are allowed.
  template<class OptData>
  void addEqualityNonLinearConstraint(OptDataConstraintFn<OptData> constraint_function,
                                      OptData data,
                                      double tolerance) {

    auto stored = std::make_shared<OptData>(std::move(data));
    ConstraintFunction closure = [constraint_function = std::move(constraint_function), stored](const std::vector<double>& x, std::vector<double>& grad) -> double {
      if (not grad.empty()) {

        ExecEnv::log().error("Optimize::addEqualityNonLinearConstraint<OptData>; expected non-differential algorithm");
        return std::numeric_limits<double>::quiet_NaN();

      }
      return constraint_function(x, *stored);
    };
    equality_constraints_.push_back({std::move(closure), tolerance});

  }

  // Convenience overloads for the common std::vector<double> data case, so that plain function pointers
  // (and other callables) can be passed without explicitly wrapping them in std::function.
  void addEqualityNonLinearConstraint(OptDerivDataConstraintFn<std::vector<double>> constraint_function,
                                      const std::vector<double>& data,
                                      double tolerance) {
    addEqualityNonLinearConstraint<std::vector<double>>(std::move(constraint_function), data, tolerance);
  }
  void addEqualityNonLinearConstraint(OptDataConstraintFn<std::vector<double>> constraint_function,
                                      const std::vector<double>& data,
                                      double tolerance) {
    addEqualityNonLinearConstraint<std::vector<double>>(std::move(constraint_function), data, tolerance);
  }

  void addEqualityNonLinearConstraint(OptDerivConstraintFn constraint_function, double tolerance);
  void addEqualityNonLinearConstraint(OptConstraintFn constraint_function, double tolerance);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // The non-linear inequality functions.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // With derivative info and data. The data is *copied* into the closure so temporaries are allowed.
  template<class OptData>
  void addInequalityNonLinearConstraint(OptDerivDataConstraintFn<OptData> constraint_function,
                                         OptData data,
                                         double tolerance) {

    auto stored = std::make_shared<OptData>(std::move(data));
    ConstraintFunction closure = [constraint_function = std::move(constraint_function), stored](const std::vector<double>& x, std::vector<double>& grad) -> double {
      return constraint_function(x, grad, *stored);
    };
    inequality_constraints_.push_back({std::move(closure), tolerance});

  }

  // Without derivative, with data. The data is *copied* into the closure so temporaries are allowed.
  template<class OptData>
  void addInequalityNonLinearConstraint(OptDataConstraintFn<OptData> constraint_function,
                                         OptData data,
                                         double tolerance) {

    auto stored = std::make_shared<OptData>(std::move(data));
    ConstraintFunction closure = [constraint_function = std::move(constraint_function), stored](const std::vector<double>& x, std::vector<double>& grad) -> double {
      if (not grad.empty()) {

        ExecEnv::log().error("Optimize::addInequalityNonLinearConstraint<OptData>; expected non-differential algorithm");
        return std::numeric_limits<double>::quiet_NaN();

      }
      return constraint_function(x, *stored);
    };
    inequality_constraints_.push_back({std::move(closure), tolerance});

  }

  // Convenience overloads for the common std::vector<double> data case, so that plain function pointers
  // (and other callables) can be passed without explicitly wrapping them in std::function.
  void addInequalityNonLinearConstraint(OptDerivDataConstraintFn<std::vector<double>> constraint_function,
                                        const std::vector<double>& data,
                                        double tolerance) {
    addInequalityNonLinearConstraint<std::vector<double>>(std::move(constraint_function), data, tolerance);
  }
  void addInequalityNonLinearConstraint(OptDataConstraintFn<std::vector<double>> constraint_function,
                                        const std::vector<double>& data,
                                        double tolerance) {
    addInequalityNonLinearConstraint<std::vector<double>>(std::move(constraint_function), data, tolerance);
  }

  void addInequalityNonLinearConstraint(OptDerivConstraintFn constraint_function, double tolerance);
  void addInequalityNonLinearConstraint(OptConstraintFn constraint_function, double tolerance);

  // Temp test function (to be removed).
  static void opt_test();

private:

  // Optimization parameters.
  OptimizationAlgorithm opt_alg_;
  std::size_t dimension_;
  OptimizationType opt_type_;
  ObjectiveFunction objective_;
  std::vector<double> upper_bound_;
  std::vector<double> lower_bound_;
  std::vector<NonLinearConstraint> equality_constraints_;
  std::vector<NonLinearConstraint> inequality_constraints_;
  std::vector<OptimalStopping> stopping_vector_;

  // Private entry to the underlying optimizer code. Stores the objective closure on this object.
  [[nodiscard]] OptResultTuple runOptimize(std::vector<double>& parameter_x_vector, ObjectiveFunction objective);

  // The returned integral type is cast to an nlopt:: optimization algorithm enum.
  [[nodiscard]] static std::size_t convertAlgorithm(OptimizationAlgorithm algorithm);

  // NLopt C-API trampolines: receive void* pointing at this Optimize (objective) or at a NonLinearConstraint.
  [[nodiscard]] static double objectiveCallback(const std::vector<double>& x, std::vector<double>& grad, void* void_self_ptr) {

    auto* self = static_cast<Optimize*>(void_self_ptr);
    return self->objective_(x, grad);

  }

  [[nodiscard]] static double constraintCallback(const std::vector<double>& x, std::vector<double>& grad, void* void_constraint_ptr) {

    auto* constraint = static_cast<NonLinearConstraint*>(void_constraint_ptr);
    return constraint->constraint_function(x, grad);

  }

};


} // namespace kellerberrin

#endif // KEL_OPTIMIZE_H