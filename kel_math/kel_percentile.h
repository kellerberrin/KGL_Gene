//
// Created by kellerberrin on 22/6/20.
//

#ifndef KEL_PERCENTILE_H
#define KEL_PERCENTILE_H

#include "kel_exec_env.h"

#include <map>
#include <vector>
#include <optional>
#include <algorithm>
#include <cmath>

namespace kellerberrin {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A template object that allows percentile stratification of payload objects based on a sortable key.
// Which will generally be an integer or double (or anything else that can be sorted).
// Currently implemented as a sorted vector.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Sortable, class Payload>
using ComparePercentile = std::function<bool(const std::pair<Sortable, Payload> &a, const std::pair<Sortable, Payload> &b)>;

template <typename Sortable, class Payload>
class Percentile {

public:

  explicit Percentile(ComparePercentile<Sortable, Payload> compare) : compare_(compare) {}
  Percentile() : compare_(default_compare_) {}
  ~Percentile() = default;

  // Adds an element to the sortable vector and sets the sort needed flag.
  // The vector is only sorted when accessed by one of the functions below.
  void addElement(Sortable sortable, Payload payload)
    { percentile_vector_.emplace_back(std::move(sortable), std::move(payload)); need_sort_ = true; }

  // Get the vector value that corresponds to the percentile. std::nullopt if an empty vector.
  [[nodiscard]] std::optional<std::pair<Sortable, Payload>> percentile(double percentile_value) const;

  // Get a range of the sorted vector that corresponds to the percentiles, i.e. (0.0, 1.0) returns the entire vector.
  [[nodiscard]] std::vector<std::pair<Sortable, Payload>> getPercentileRange(double lower_percentile, double upper_percentile) const;

  // ReturnType the sorted percentile vector.
  [[nodiscard]] const std::vector<std::pair<Sortable, Payload>>& getVector() const { conditionalSort(); return percentile_vector_; }

  // Given a sortable value find the corresponding percentile, 0.0 if an empty vector. Payload is not used.
  [[nodiscard]] double inversePercentile(Sortable value, Payload dummy_payload) const;

  // Given a value, find the number of elements in the array >= to the value.
  [[nodiscard]] size_t findGEQCount(Sortable find_value, Payload dummy_payload) const;

private:

  // A const access may require the vector to be re-sorted (if it has been updated).
  mutable std::vector<std::pair<Sortable, Payload>> percentile_vector_;
  mutable bool need_sort_{true};

  // Default sort lambda, generally used at runtime, although users can specify a bespoke sort function.
  constexpr static auto default_compare_
     = []( const std::pair<Sortable, Payload> &a,
           const std::pair<Sortable, Payload> &b) -> bool { return a.first < b.first; };

  // Runtime sort function, see above.
  ComparePercentile<Sortable, Payload> compare_;

  // ReturnType the vector element corresponding to a percentile.
  [[nodiscard]] size_t index(double percentile) const;
  // Sort the vector if it has been updated.
  void conditionalSort() const;
  // Given a value, find the array index.
  [[nodiscard]] std::optional<size_t> findArrayIndex(Sortable find_value, Payload dummy_payload) const;

};

// ReturnType a vector range between two percentiles.
template <typename Sortable, class Payload>
std::vector<std::pair<Sortable, Payload>> Percentile<Sortable, Payload>::getPercentileRange(double lower_percentile, double upper_percentile) const {

  if (percentile_vector_.empty()) {

    return std::vector<std::pair<Sortable, Payload>>();

  }

  conditionalSort();

  auto lower_range = percentile_vector_.begin() + index(lower_percentile);
  auto upper_range = percentile_vector_.begin() + index(upper_percentile);

  return std::vector<std::pair<Sortable, Payload>>(lower_range, upper_range);

}

// ReturnType the record corresponding to the percentile.
template <typename Sortable, class Payload>
std::optional<std::pair<Sortable, Payload>> Percentile<Sortable, Payload>::percentile(double percentile_value) const {

  if (percentile_vector_.empty()) {

    return std::nullopt;

  }

  conditionalSort();

  auto iterator = percentile_vector_.begin() + index(percentile_value);

  return *iterator;

}

// Given a sortable value find the corresponding percentile, 0.0 if an empty vector. Payload is not used.
template <typename Sortable, class Payload>
double Percentile<Sortable, Payload>::inversePercentile(Sortable value, Payload dummy_payload) const {

  std::optional<size_t> index_opt = findArrayIndex(value, dummy_payload);

  if (not index_opt) {

    return 0.0;

  }

  if (percentile_vector_.size() == 1) {

    return 1.0;

  }

  return static_cast<double>(index_opt.value()) / (static_cast<double>(percentile_vector_.size()) - 1.0);

}

template <typename Sortable, class Payload>
size_t Percentile<Sortable, Payload>::findGEQCount(Sortable find_value, Payload dummy_payload) const {

  std::optional<size_t> index_opt = findArrayIndex(find_value, dummy_payload);

  if (not index_opt) {

    return 0;

  }

  return percentile_vector_.size() - index_opt.value();

}


// Calculate the vector index for a given percentile.
template <typename Sortable, class Payload>
size_t Percentile<Sortable, Payload>::index(double percentile) const {

  if (percentile < 0 or percentile > 1) {

    ExecEnv::log().error("Percentile::index, specified percentile value: {} is out of range", percentile);
    return 0;

  }

  if (percentile_vector_.empty()) {

    ExecEnv::log().error("Percentile::index, attempt to generate an index of an empty vector");
    return 0;

  }

  // Calculate the vector index.
  double proportion = (static_cast<double>(percentile_vector_.size()) * percentile) - 0.5;
  long vector_index = std::lround(proportion);
  vector_index = vector_index < 0 ? 0 : vector_index;
  vector_index = vector_index > static_cast<long>(percentile_vector_.size() - 1) ? (percentile_vector_.size() - 1) : vector_index;

  return static_cast<size_t>(vector_index);

}

// Binary lookup on the array index.
template <typename Sortable, class Payload>
std::optional<size_t> Percentile<Sortable, Payload>::findArrayIndex(Sortable find_value, Payload dummy_payload) const {

  if (percentile_vector_.empty()) {

    return std::nullopt;

  }

  conditionalSort();

  auto result = std::lower_bound( percentile_vector_.begin(),
                                  percentile_vector_.end(),
                                  std::pair<Sortable, Payload>(find_value, dummy_payload),
                                  compare_);

  if (result == percentile_vector_.end()) {

    return percentile_vector_.size() -1; // last element

  }

  return std::distance(percentile_vector_.begin(), result);

}


// If Updated, sort the percentile vector before any access.
template <typename Sortable, class Payload>
void Percentile<Sortable, Payload>::conditionalSort() const {

  if (need_sort_) {

    std::sort(percentile_vector_.begin(), percentile_vector_.end(), compare_);

    need_sort_ = false;

  }

}


} // namespace


#endif // KEL_PERCENTILE_H
