//
// Created by kellerberrin on 18/9/20.
//

#ifndef KEL_OPTIMIZE_H
#define KEL_OPTIMIZE_H

#include <vector>
#include <string>
#include <functional>

namespace kellerberrin {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The optimize object is a facade in front of the nlopt optimization library
// If this object is used then the executable must link to "libnlopt".
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Optimization Algorithms. A mix of constrained local and global algorithms, some require
// derivatives to speed convergence. See the nlopt documentation for further description.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  LN_NELDERMEAD,
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

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Optimization return results. Same values as returned by the nlopt package, see documentation.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// How the objective function is optimized
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class OptimizationType {
  MAXIMIZE = 1,
  MINIMIZE = 2
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Stopping criteria (can be more than 1 criteria)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class OptimizeStoppingType {
  FUNCTION_VALUE = 1,     // Stop when a particular function value is reached.
  RELATIVE_FUNCTION_THRESHOLD = 2,    // Stop when the relative objective function update is below a threshold.
  ABSOLUTE_FUNCTION_THRESHOLD = 3,    // Stop when the absolute objective function update is below a threshold.
  RELATIVE_PARAMETER_THRESHOLD = 4,  // Stop when the relative weighted (normed) parameter vector update is below a threshold.
  RELATIVE_PARAMETER_WEIGHTS = 5,  // Set the relative parameter norm weights.
  ABSOLUTE_PARAMETER_THRESHOLD = 6,  // Stop when the parameter vector update is below a threshold.
  MAXIMUM_EVALUATIONS = 7, // Stop when the maximum number of evaluations have been reached.
  MAXIMUM_TIME = 8    // Stop when the specified time in seconds has elapsed.
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The optimize class, this a thin facade over the corresponding nlopt functionality.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ObjectiveFunction =  std::function<double(const std::vector<double> &x, std::vector<double> &grad, void* f_data)>;
using NonLinearConstraintFunction = ObjectiveFunction;

class Optimize {

public:


  Optimize( OptimizationAlgorithm opt_alg,
            size_t dimension,
            OptimizationType opt_type,
            ObjectiveFunction objective)
    : opt_alg_(opt_alg), dimension_(dimension), opt_type_(opt_type),  objective_(std::move(objective)) {}
  ~Optimize() = default;

  // Perform the optimization. The initial parameter vector is updated to the function x parameter stopping value.
  // The returned tuple is [result_code, optimal_function_value, iterations]
  std::tuple<OptimizationResult, double, size_t> optimize(std::vector<double>& x_parameter_vector);
  // Convert the return value to a string
  static std::string returnDescription(OptimizationResult result);
  // Define the hypercube which contains the solution. Must be the same dimension as the objective function.
  // Empty vector is ignored
  void boundingHypercube(const std::vector<double>& upper_bound = {}, const std::vector<double>& lower_bound = {});
  // Equality Constraints
  void addEqualityNonLinearConstraint(NonLinearConstraintFunction constraint_function, const std::vector<double>& data, double tolerance);
  // Inequality Constraints
  void addInequalityNonLinearConstraint(NonLinearConstraintFunction constraint_function, const std::vector<double>& data, double tolerance);
  // Stopping Criteria. The stopping criteria vector either has 1 element.
  // Or for x parameter stopping criteria (only), the same dimensionality of the objective function parameters.
  void stoppingCriteria(OptimizeStoppingType stopping_type, const std::vector<double>& stopping_value);

  static void opt_test();

private:

  OptimizationAlgorithm opt_alg_;
  size_t dimension_;
  OptimizationType opt_type_;
  ObjectiveFunction objective_;
  std::vector<double> upper_bound_;
  std::vector<double> lower_bound_;
  struct NonLinearConstraint {
    NonLinearConstraintFunction constraint_function;
    std::vector<double> data;
    double tolerance;
  };
  std::vector<NonLinearConstraint> equality_constraints_;
  std::vector<NonLinearConstraint> inequality_constraints_;
  struct OptimalStopping {
    OptimizeStoppingType stopping_type;
    const std::vector<double> stopping_value;
  };
  std::vector<OptimalStopping> stopping_vector_;

  // The returned integral type is cast to an nlopt:: optimization algorithm enum.
  static size_t convertAlgorithm(OptimizationAlgorithm algorithm);

};


} // namespace

#endif //KEL_OPTIMIZE_H
