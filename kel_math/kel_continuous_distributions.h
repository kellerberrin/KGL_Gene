// Kellerberrin 2026.

#ifndef KEL_CONTINUOUS_DISTRIBUTIONS_H
#define KEL_CONTINUOUS_DISTRIBUTIONS_H

#include <algorithm>
#include <random>

#include "kel_entropy_source.h"
#include "kel_exec_env.h"

namespace kellerberrin {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Continuous distributions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Real (double) random number generator on the interval [lower_bound, upper_bound]
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UniformRealDistribution final {
public:

  explicit UniformRealDistribution(double lower_bound, double upper_bound) {
    if (lower_bound > upper_bound) {

      ExecEnv::log().warn("UniformRealDistribution; lower_bound:{} exceeds upper_bound:{}, swapping bounds",
                          lower_bound, upper_bound);
      std::swap(lower_bound, upper_bound);

    }
    dist_ = std::uniform_real_distribution<>(lower_bound, upper_bound);
  }

  UniformRealDistribution(const UniformRealDistribution&) = delete;
  UniformRealDistribution(UniformRealDistribution&&) = default;
  UniformRealDistribution& operator=(const UniformRealDistribution&) = delete;
  UniformRealDistribution& operator=(UniformRealDistribution&&) = default;
  ~UniformRealDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] double random(URBG& source) const { return dist_(source); }

private:

  mutable std::uniform_real_distribution<> dist_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Random numbers on the unit interval [0, 1]
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UniformUnitDistribution final {
public:

  UniformUnitDistribution() = default;

  UniformUnitDistribution(const UniformUnitDistribution&) = delete;
  UniformUnitDistribution(UniformUnitDistribution&&) = default;
  UniformUnitDistribution& operator=(const UniformUnitDistribution&) = delete;
  UniformUnitDistribution& operator=(UniformUnitDistribution&&) = default;
  ~UniformUnitDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] double random(URBG& source) const { return dist_(source); }

private:

  mutable std::uniform_real_distribution<> dist_{0.0, 1.0};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Normally distributed random numbers with mean and std_deviation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NormalDistribution final {
public:

  explicit NormalDistribution(double mean, double std_deviation)
        : dist_(mean, std_deviation) {}

  NormalDistribution(const NormalDistribution&) = delete;
  NormalDistribution(NormalDistribution&&) = default;
  NormalDistribution& operator=(const NormalDistribution&) = delete;
  NormalDistribution& operator=(NormalDistribution&&) = default;
  ~NormalDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] double random(URBG& source) const { return dist_(source); }

  [[nodiscard]] static double pdf(double x, double mean, double std_deviation);

private:

  mutable std::normal_distribution<> dist_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard normally distributed random numbers with mean = 0.0 and std_deviation = 1.0.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class StdNormalDistribution final {
public:

  StdNormalDistribution() = default;
  StdNormalDistribution(const StdNormalDistribution&) = delete;
  StdNormalDistribution(StdNormalDistribution&&) = default;
  StdNormalDistribution& operator=(const StdNormalDistribution&) = delete;
  StdNormalDistribution& operator=(StdNormalDistribution&&) = default;
  ~StdNormalDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] double random(URBG& source) const { return dist_(source); }

private:

  mutable std::normal_distribution<> dist_{0.0, 1.0};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gamma distribution.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GammaDistribution final {
public:

  explicit GammaDistribution(double shape, double scale)
        : dist_(shape, scale), shape_(shape), scale_(scale) {}

  GammaDistribution(const GammaDistribution&) = delete;
  GammaDistribution(GammaDistribution&&) = default;
  GammaDistribution& operator=(const GammaDistribution&) = delete;
  GammaDistribution& operator=(GammaDistribution&&) = default;
  ~GammaDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] double random(URBG& source) const { return dist_(source); }

  [[nodiscard]] double pdf(double x) const;
  [[nodiscard]] double cdf(double x) const;
  [[nodiscard]] double quantile(double p) const;

  [[nodiscard]] double shape() const noexcept { return shape_; }
  [[nodiscard]] double scale() const noexcept { return scale_; }

private:

  mutable std::gamma_distribution<> dist_;
  double shape_;
  double scale_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta distributed random numbers implemented as a ratio of gamma random variates.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BetaDistribution final {
public:

  explicit BetaDistribution(double a, double b)
        : x_(a, 1.0), y_(b, 1.0), a_(a), b_(b) {}

  BetaDistribution(const BetaDistribution&) = delete;
  BetaDistribution(BetaDistribution&&) = default;
  BetaDistribution& operator=(const BetaDistribution&) = delete;
  BetaDistribution& operator=(BetaDistribution&&) = default;
  ~BetaDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] double random(URBG& source) const {
    const double x = x_(source);
    const double y = y_(source);
    return x / (x + y);
  }

  [[nodiscard]] double a() const noexcept { return a_; }
  [[nodiscard]] double b() const noexcept { return b_; }

  [[nodiscard]] static double logInverseBetaFunction(double a, double b);
  [[nodiscard]] static double logPartialPdf(double x, double a, double b);
  [[nodiscard]] static double logPdf(double x, double a, double b);
  [[nodiscard]] static double pdf(double x, double a, double b);
  [[nodiscard]] static double mean(double a, double b);
  [[nodiscard]] static double var(double a, double b);
  [[nodiscard]] static double mode(double a, double b);

private:

  mutable std::gamma_distribution<> x_;
  mutable std::gamma_distribution<> y_;
  double a_;
  double b_;
};


}   // namespace kellerberrin

#endif // KEL_CONTINUOUS_DISTRIBUTIONS_H