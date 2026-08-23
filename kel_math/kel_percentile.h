// Kellerberrin 2026.

#ifndef KEL_PERCENTILE_H
#define KEL_PERCENTILE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include "kel_exec_env.h"

namespace kellerberrin {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Percentile stratification of (key, payload) pairs.
//
// Sortable: type with a strict weak ordering (e.g. int, double).
// Payload:  arbitrary attached value.
// Compare:  callable bool(const Sortable&, const Sortable&).  Defaults to std::less<>.
//
// The container is lazily sorted: addElement() only marks data dirty; actual sorting happens on first read.
//
// IMPORTANT: This class is NOT thread-safe.  Any reference returned by getVector() is invalidated by addElement().
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Sortable, typename Payload, typename Compare = std::less<>>
class Percentile {
public:

  using value_type = std::pair<Sortable, Payload>;

  explicit Percentile(Compare compare = Compare{}) : compare_(std::move(compare)) {}

  ~Percentile() = default;

  Percentile(const Percentile&) = default;
  Percentile& operator=(const Percentile&) = default;
  Percentile(Percentile&&) = default;
  Percentile& operator=(Percentile&&) = default;

  // Add a (key, payload) pair.  Invalidates any reference previously returned by getVector().
  void addElement(Sortable sortable, Payload payload) {
    data_.emplace_back(std::move(sortable), std::move(payload));
    need_sort_ = true;
  }

  // Return the element closest to the requested percentile, or std::nullopt if empty.
  // Percentile must be in [0.0, 1.0]; out-of-range and empty cases are logged and return std::nullopt.
  [[nodiscard]] std::optional<value_type> percentile(double percentile_value) const;

  // Return the inclusive range [lower_percentile, upper_percentile].
  // Returns an empty vector if the container is empty or if lower > upper.
  [[nodiscard]] std::vector<value_type> getPercentileRange(double lower_percentile,
                                                            double upper_percentile) const;

  // Return a const reference to the sorted underlying vector.
  [[nodiscard]] const std::vector<value_type>& getVector() const;

  // Return the percentile rank of a key value.  Returns 0.0 if empty.
  // For a single-element container, returns 1.0 (matches legacy behavior).
  [[nodiscard]] double inversePercentile(const Sortable& value) const;

  // Return the number of elements with key >= value.
  [[nodiscard]] std::size_t findGEQCount(const Sortable& value) const;

  [[nodiscard]] bool empty() const noexcept { return data_.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }

  void clear() noexcept {
    data_.clear();
    need_sort_ = false;
  }

private:

  mutable std::vector<value_type> data_;
  mutable bool need_sort_{true};
  Compare compare_;

  [[nodiscard]] std::size_t index(double percentile) const;
  void ensureSorted() const;
};


template <typename Sortable, typename Payload, typename Compare>
void Percentile<Sortable, Payload, Compare>::ensureSorted() const {

  if (need_sort_) {

    std::sort(data_.begin(), data_.end(),
              [this](const value_type& a, const value_type& b) {
                return compare_(a.first, b.first);
              });
    need_sort_ = false;

  }

}

template <typename Sortable, typename Payload, typename Compare>
std::size_t Percentile<Sortable, Payload, Compare>::index(double percentile) const {

  if (percentile < 0.0 || percentile > 1.0) {

    ExecEnv::log().error("Percentile value: {} must be in [0.0, 1.0]", percentile);
    return 0;

  }

  if (data_.empty()) {

    ExecEnv::log().warn("Cannot compute percentile index of an empty distribution");
    return 0;

  }

  // Preserve the original index formula: round(N * p - 0.5), clamped.
  const double proportion = (static_cast<double>(data_.size()) * percentile) - 0.5;
  long long vector_index = std::llround(proportion);
  vector_index = std::clamp(vector_index,
                            0LL,
                            static_cast<long long>(data_.size() - 1));

  return static_cast<std::size_t>(vector_index);

}

template <typename Sortable, typename Payload, typename Compare>
std::optional<std::pair<Sortable, Payload>>
Percentile<Sortable, Payload, Compare>::percentile(double percentile_value) const {

  if (data_.empty()) {

    return std::nullopt;

  }

  ensureSorted();
  return data_[index(percentile_value)];

}

template <typename Sortable, typename Payload, typename Compare>
std::vector<std::pair<Sortable, Payload>>
Percentile<Sortable, Payload, Compare>::getPercentileRange(double lower_percentile,
                                                            double upper_percentile) const {

  if (lower_percentile < 0.0 || upper_percentile > 1.0) {

    ExecEnv::log().error("Percentile interval: [{}, {}] must be in [0.0, 1.0]",
                         lower_percentile, upper_percentile);
    return {};

  }

  if (data_.empty() || lower_percentile > upper_percentile) {

    return {};

  }

  ensureSorted();

  const std::size_t lower_idx = index(lower_percentile);
  const std::size_t upper_idx = std::min(index(upper_percentile) + 1, data_.size());

  return std::vector<value_type>(data_.begin() + lower_idx, data_.begin() + upper_idx);

}

template <typename Sortable, typename Payload, typename Compare>
const std::vector<std::pair<Sortable, Payload>>&
Percentile<Sortable, Payload, Compare>::getVector() const {

  ensureSorted();
  return data_;

}

template <typename Sortable, typename Payload, typename Compare>
double Percentile<Sortable, Payload, Compare>::inversePercentile(const Sortable& value) const {

  if (data_.empty()) {

    return 0.0;

  }

  ensureSorted();

  auto it = std::lower_bound(data_.begin(), data_.end(), value,
                              [this](const value_type& pair, const Sortable& key) {
                                return compare_(pair.first, key);
                              });

  if (data_.size() == 1) {

    return 1.0;

  }

  const std::size_t idx = std::min(static_cast<std::size_t>(std::distance(data_.begin(), it)),
                                   data_.size() - 1);
  return static_cast<double>(idx) / static_cast<double>(data_.size() - 1);

}

template <typename Sortable, typename Payload, typename Compare>
std::size_t Percentile<Sortable, Payload, Compare>::findGEQCount(const Sortable& value) const {

  if (data_.empty()) {

    return 0;

  }

  ensureSorted();

  auto it = std::lower_bound(data_.begin(), data_.end(), value,
                              [this](const value_type& pair, const Sortable& key) {
                                return compare_(pair.first, key);
                              });

  return static_cast<std::size_t>(std::distance(it, data_.end()));

}

} // namespace kellerberrin

#endif // KEL_PERCENTILE_H