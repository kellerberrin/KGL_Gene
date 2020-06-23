//
// Created by kellerberrin on 22/6/20.
//

#ifndef KEL_PERCENTILE_H
#define KEL_PERCENTILE_H

#include "kel_exec_env.h"

#include <map>
#include <vector>
#include <algorithm>

namespace kellerberrin {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A template object that allows percentile stratification of payload objects based on a sortable key.
// Which will generally be an integer or double (or anything else that can be sorted).
// Currently implemented as a sorted vector.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Sortable, class Payload>
class Percentile {

public:

  Percentile() = default;
  ~Percentile() = default;

  void addElement(Sortable sortable, Payload payload)
    { percentile_vector_.emplace_back(std::move(sortable), std::move(payload)); need_sort_ = true; }

  std::optional<std::pair<Sortable, Payload>> percentile(double percentile_value) const;

  std::vector<std::pair<Sortable, Payload>> getPercentileRange(double lower_percentile, double upper_percentile) const;

  const std::vector<std::pair<Sortable, Payload>>& getVector() const { conditionalSort(); return percentile_vector_; }

private:

  mutable std::vector<std::pair<Sortable, Payload>> percentile_vector_;
  mutable bool need_sort_{true};

  [[nodiscard]] size_t index(double percentile) const;
  void conditionalSort() const;

};

// Return a vector range between two percentiles.
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

// Return the record corresponding to the percentile.
template <typename Sortable, class Payload>
std::optional<std::pair<Sortable, Payload>> Percentile<Sortable, Payload>::percentile(double percentile_value) const {

  if (percentile_vector_.empty()) {

    return std::nullopt;

  }

  conditionalSort();

  auto iterator = percentile_vector_.begin() + index(percentile_value);

  return *iterator;

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

// If Updated, sort the percentile vector before any access.
template <typename Sortable, class Payload>
void Percentile<Sortable, Payload>::conditionalSort() const {

  if (need_sort_) {

    auto comp = [](const std::pair<Sortable, Payload> &a, const std::pair<Sortable, Payload> &b) -> bool {
      return a.first < b.first;
    };

    std::sort(percentile_vector_.begin(), percentile_vector_.end(), comp);

    need_sort_ = false;

  }

}



} // namespace


#endif // KEL_PERCENTILE_H
