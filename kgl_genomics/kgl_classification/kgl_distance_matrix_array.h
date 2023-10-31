//
// Created by kellerberrin on 31/10/23.
//

#ifndef KGL_DISTANCE_MATRIX_ARRAY_H
#define KGL_DISTANCE_MATRIX_ARRAY_H



#include "kel_exec_env.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of a strict lower triangular calculateDistance matrix using the Boost library.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization level namespace


template <typename T>
concept NumericTypeArray = std::integral<T> || std::floating_point<T>;

template <NumericTypeArray MatrixValue>
class DistanceMatrixArray {

  [[nodiscard]] static size_t lower_triangle_size(size_t n) { return ((n+1) * n) / 2; };
  [[nodiscard]] size_t lower_triangle_index(size_t row, size_t column) const {

    if (column > row) std::swap(row, column);

    return row * matrix_size_ + ((row * (row-1)) / 2) + column;

  };

  struct AuxMatrixArray{

    size_t row_;
    size_t column_;
    MatrixValue value_;

  };

  std::vector<AuxMatrixArray> createAuxArray() const {

    std::vector<AuxMatrixArray> aux_array(lower_triangular_.size());

    for (size_t idx = 0; idx < lower_triangular_.size(); ++idx) {

      aux_array[idx].value_ = lower_triangular_[idx];

    }

    for (size_t row = 0; row < matrix_size_; ++row) {

      for (size_t column = row; column < matrix_size_; ++column) {

        size_t vector_index = lower_triangle_index(row, column);
        aux_array[vector_index].row_ = row;
        aux_array[vector_index].column_ = column;

      }

    }

    return aux_array;

  }

public:

  explicit DistanceMatrixArray(size_t matrix_size)  { resize(matrix_size); }
  DistanceMatrixArray(const DistanceMatrixArray &) = default;
  ~DistanceMatrixArray() = default;

  [[nodiscard]] size_t size() const { return matrix_size_; }
  void resize(size_t new_size) {

    size_t vector_size = lower_triangle_size(new_size);
    lower_triangular_.resize( vector_size, 0.0 );
    matrix_size_ = new_size;

  }

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

  std::vector<MatrixValue> lower_triangular_;
  size_t matrix_size_;

};

template <NumericTypeArray MatrixValue>
MatrixValue DistanceMatrixArray<MatrixValue>::getDistance(size_t row, size_t column) const {

  if (row >= size() or column >= size()) {

    ExecEnv::log().error("Index too large row: {}, column:{} for calculateDistance matrix size: {}", row, column, size());
    return 0;

  }

  size_t vector_index = lower_triangle_index(row, column);

  if (vector_index >= lower_triangular_.size()) {

    ExecEnv::log().error("Index too large: {}, for vector size: {}", vector_index, lower_triangular_.size());
    return 0;

  }

  return lower_triangular_[vector_index];

}

template <NumericTypeArray MatrixValue>
void DistanceMatrixArray<MatrixValue>::setDistance(size_t row, size_t column, MatrixValue distance) {

  if (row >= size() or column >= size()) {

    ExecEnv::log().error("Index too large row: {}, column:{} for calculateDistance matrix size: {}", row, column, size());
    return;

  }

  size_t vector_index = lower_triangle_index(row, column);

  if (vector_index >= lower_triangular_.size()) {

    ExecEnv::log().error("Index too large: {}, for vector size: {}", vector_index, lower_triangular_.size());
    return;

  }

  lower_triangular_[vector_index] = distance;

}

// Sort search with a cost of O(n * log n)
template <NumericTypeArray MatrixValue>
std::tuple<MatrixValue, size_t, size_t> DistanceMatrixArray<MatrixValue>::minimum() const {

  if (lower_triangular_.empty()) {

    ExecEnv::log().error("Matrix vector is zero sized");
    return { 0.0, 0, 0};

  }

  auto aux_array = createAuxArray();
  auto sort_comp = [] ( const AuxMatrixArray& lhs, const AuxMatrixArray& rhs ) { return lhs.value_ < rhs.value_; };
  std::sort(aux_array.begin(), aux_array.end(), sort_comp);
  auto& min = aux_array.front();

  return {min.value_, min.row_, min.column_};

}

// Sort search with a cost of O(n * log n)
template <NumericTypeArray MatrixValue>
std::tuple<MatrixValue, size_t, size_t> DistanceMatrixArray<MatrixValue>::maximum() const {

  if (lower_triangular_.empty()) {

    ExecEnv::log().error("Matrix vector is zero sized");
    return { 0.0, 0, 0};

  }

  auto aux_array = createAuxArray();
  auto sort_comp = [] ( const AuxMatrixArray& lhs, const AuxMatrixArray& rhs ) { return lhs.value_ < rhs.value_; };
  std::sort(aux_array.begin(), aux_array.end(), sort_comp);
  auto& min = aux_array.back();

  return {min.value_, min.row_, min.column_};

}


// Sort search with a cost of O(n * log n)
template <NumericTypeArray MatrixValue>
std::pair<MatrixValue, MatrixValue> DistanceMatrixArray<MatrixValue>::max_min() const {

  std::vector<MatrixValue> vector_copy(lower_triangular_);
  std::sort(vector_copy.begin(), vector_copy.begin());

  if (vector_copy.empty()) {

    ExecEnv::log().error("Matrix vector is zero sized");
    return { 0.0, 0.0};

  }

  return {vector_copy.back(), vector_copy.front()};

}

template <NumericTypeArray MatrixValue>
std::pair<MatrixValue, MatrixValue> DistanceMatrixArray<MatrixValue>::max_min_diagonal() const {

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

template <NumericTypeArray MatrixValue>
void DistanceMatrixArray<MatrixValue>::normalizeDistance() {

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
template <NumericTypeArray MatrixValue>
void DistanceMatrixArray<MatrixValue>::zeroDiagonal() {

  auto [max, min] = max_min_diagonal();
  if (max != 0.0 or min != 0.0) {

    ExecEnv::log().error("Prior matrix diagonal range: [{}, {}]; diagonals now set to 0.0", max, min);

    for (size_t row = 0; row < size(); ++row) {

      setDistance(row, row, 0.0);

    }

  }

}

} // Namespace

#endif //KGL_DISTANCE_MATRIX_ARRAY_H
