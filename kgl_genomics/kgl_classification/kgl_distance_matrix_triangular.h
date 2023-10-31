//
// Created by kellerberrin on 31/10/23.
//

#ifndef KGL_DISTANCE_MATRIX_TRIANGULAR_H
#define KGL_DISTANCE_MATRIX_TRIANGULAR_H


#include "kel_exec_env.h"

#include <boost/numeric/ublas/triangular.hpp>

namespace bnu = boost::numeric::ublas;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of a strict lower triangular calculateDistance matrix using the Boost library.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization level namespace


template <typename T>
concept NumericTypeTriangle = std::integral<T> || std::floating_point<T>;

template <NumericTypeTriangle MatrixValue>
class DistanceMatrixTriangular {

  using DistanceMatrixImplType = bnu::triangular_matrix<MatrixValue, bnu::strict_lower>;

public:

  explicit DistanceMatrixTriangular(size_t matrix_size) : lower_triangular_(matrix_size, matrix_size) {}

  explicit DistanceMatrixTriangular(const DistanceMatrixTriangular &) = default;
  ~DistanceMatrixTriangular() = default;

  [[nodiscard]] size_t size() const { return lower_triangular_.size1(); }
  void resize(size_t new_size) { lower_triangular_.resize(new_size, new_size, false); }

  [[nodiscard]] MatrixValue getDistance(size_t i, size_t j) const;
  void setDistance(size_t i, size_t j, MatrixValue distance);

  // Tuple returns value, row index, column index in that order.
  [[nodiscard]] std::tuple<MatrixValue, size_t, size_t> minimum() const;
  [[nodiscard]] std::tuple<MatrixValue, size_t, size_t> maximum() const;

  // the .first element is the maximum value, the .second element is the minmum value.
  [[nodiscard]] std::pair<MatrixValue, MatrixValue> max_min() const;
  // the .first element is the maximum value, the .second element is the minmum value.
  [[nodiscard]] std::pair<MatrixValue, MatrixValue> max_min_diagonal() const;

  // Rescale elements to the interval [0, 1].
  void normalizeDistance();
  // Zero all diagonal elements.
  void zeroDiagonal();

private:

  DistanceMatrixImplType lower_triangular_;

};

template <NumericTypeTriangle MatrixValue>
MatrixValue DistanceMatrixTriangular<MatrixValue>::getDistance(size_t i, size_t j) const {

  if (i >= size() or i >= size()) {

    ExecEnv::log().error("Index too large i: {}, j:{} for calculateDistance matrix size: {}", i, j, size());
    return 0;

  }

  if (j > i) {

    std::swap(i, j);

  }

  return lower_triangular_(i, j);

}

template <NumericTypeTriangle MatrixValue>
void DistanceMatrixTriangular<MatrixValue>::setDistance(size_t i, size_t j, MatrixValue distance) {

  if (i >= size() or i >= size()) {

    ExecEnv::log().error("Index too large i: {}, j:{} for calculateDistance matrix size: {}", i, j, size());
    return;

  }

  if (j > i) {

    std::swap(i, j);

  }

  lower_triangular_(i, j) = distance;

}

// Exhaustive search with a cost of O(n^2)
template <NumericTypeTriangle MatrixValue>
std::tuple<MatrixValue, size_t, size_t> DistanceMatrixTriangular<MatrixValue>::minimum() const {

  bool first_pass{true};
  size_t i{0};
  size_t j{0};
  MatrixValue min_distance{0};
  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      if (first_pass) {

        min_distance = getDistance(row, column);
        i = row;
        j = column;
        first_pass = false;

      } else {

        auto check_min = getDistance(row, column);
        if (check_min < min_distance) {

          min_distance = check_min;
          i = row;
          j = column;

        }

      }

    }

  }

  return {min_distance, i, j};

}

// Exhaustive search with a cost of O(n^2)
template <NumericTypeTriangle MatrixValue>
std::tuple<MatrixValue, size_t, size_t> DistanceMatrixTriangular<MatrixValue>::maximum() const {

  bool first_pass = true;
  size_t i{0};
  size_t j{0};
  MatrixValue max_distance = 0;
  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      if (first_pass) {

        max_distance = getDistance(row, column);
        i = row;
        j = column;
        first_pass = false;

      } else {

        MatrixValue check_max = getDistance(row, column);
        if (check_max > max_distance) {

          max_distance = check_max;
          i = row;
          j = column;

        }

      }

    }

  }

  return {max_distance, i, j};

}


// Exhaustive search with a cost of O(n^2)
template <NumericTypeTriangle MatrixValue>
std::pair<MatrixValue, MatrixValue> DistanceMatrixTriangular<MatrixValue>::max_min() const {

  bool first_pass = true;
  MatrixValue max_distance{0};
  MatrixValue min_distance{0};
  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      if (first_pass) {

        max_distance = min_distance = getDistance(row, column);
        first_pass = false;

      } else {

        MatrixValue check_value = getDistance(row, column);
        if (check_value > max_distance) {

          max_distance = check_value;

        }
        if (check_value < min_distance) {

          min_distance = check_value;

        }

      }

    }

  }

  return {max_distance, min_distance};

}

template <NumericTypeTriangle MatrixValue>
std::pair<MatrixValue, MatrixValue> DistanceMatrixTriangular<MatrixValue>::max_min_diagonal() const {

  bool first_pass{true};
  MatrixValue max_distance{0};
  MatrixValue min_distance{0};
  for (size_t row = 0; row < size(); ++row) {

    if (first_pass) {

      max_distance = min_distance = getDistance(row, row);
      first_pass = false;

    } else {

      MatrixValue check_value = getDistance(row, row);
      if (check_value > max_distance) {

        max_distance = check_value;

      }
      if (check_value < min_distance) {

        min_distance = check_value;

      }

    }

  }

  return {max_distance, min_distance};

}

template <NumericTypeTriangle MatrixValue>
void DistanceMatrixTriangular<MatrixValue>::normalizeDistance() {

  auto [max, min] = max_min();
  MatrixValue range = max - min;

  if (range == 0.0) {

    ExecEnv::log().error("CcalculateDistance range for all nodes is zero");
    return;

  }

  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; column++) {

      MatrixValue raw_distance = getDistance(row, column);
      MatrixValue adj_distance = (raw_distance - min) / range;
      setDistance(row, column, adj_distance);

    }

  }

}

// Zero all diagonal elements.
template <NumericTypeTriangle MatrixValue>
void DistanceMatrixTriangular<MatrixValue>::zeroDiagonal() {

  auto [max, min] = max_min_diagonal();
  if (max != 0.0 or min != 0.0) {

    ExecEnv::log().error("Prior matrix diagonal range: [{}, {}]; diagonals now set to 0.0", max, min);

    for (size_t row = 0; row < size(); ++row) {

      setDistance(row, row, 0.0);

    }

  }

}

} // Namespace


#endif //KGL_DISTANCE_MATRIX_TRIANGULAR_H
