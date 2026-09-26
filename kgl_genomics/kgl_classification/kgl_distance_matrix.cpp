//
// kgl_distance_matrix.cpp — packed strict-lower triangular distance matrix.
//
// Created by kellerberrin on 30/10/23.
//


#include "kel_exec_env.h"
#include "kgl_distance_matrix.h"

#include <algorithm>
#include <ranges>
#include <vector>


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DistanceMatrixImpl — strict lower triangular matrix stored as a packed vector.
//
// The element (i, j) with i > j is stored at packed offset (i * (i - 1) / 2) + j; elements on and
// above the diagonal do not exist in storage. This replaces the previous boost::ublas triangular
// matrix (previously the only Boost dependency in the kgl_genomics layer) with identical
// observable semantics: symmetric index swapping, index and diagonal guards, and row-major
// scan order for minimum() and maximum().
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceMatrixImpl {

public:

  explicit DistanceMatrixImpl(size_t matrix_size = 1) { resize(matrix_size); }
  DistanceMatrixImpl(const DistanceMatrixImpl&) = default;
  ~DistanceMatrixImpl() = default;

  void addMatrix(const DistanceMatrixImpl& add_matrix) {

    if (add_matrix.size() != size()) {

      ExecEnv::log().error("Cannot add matrices of different sizes; lhs size: {}, rhs size: {}",
                           size(), add_matrix.size());
      return;

    }

    for (size_t idx = 0; idx < packed_matrix_.size(); ++idx) {

      packed_matrix_[idx] += add_matrix.packed_matrix_[idx];

    }

  }

  void copyMatrix(const DistanceMatrixImpl& copy_matrix) {

    matrix_size_ = copy_matrix.matrix_size_;
    packed_matrix_ = copy_matrix.packed_matrix_;

  }

  [[nodiscard]] size_t size() const { return matrix_size_; }

  void resize(size_t new_size) {

    matrix_size_ = new_size;
    packed_matrix_.assign(packedSize(new_size), 0.0);

  }

  [[nodiscard]] DistanceType_t getDistance(size_t i, size_t j) const {

    if (i >= size() or j >= size()) {

      ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
      return 0;

    }

    if (i == j) {

      ExecEnv::log().error("Matrix size: {} is strict lower triangular; [i: {}, j: {}] == 0", size(), i, j);
      return 0;

    }

    if (j > i) {

      std::swap(i, j);

    }

    return packed_matrix_[offset(i, j)];

  }

  void setDistance(size_t i, size_t j, DistanceType_t distance) {

    if (i >= size() or j >= size()) {

      ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
      return;

    }

    if (i == j) {

      // The diagonal is structural zero; there is no stored element to update.
      ExecEnv::log().error("Matrix size: {} is strict lower triangular; [i: {}, j: {}] != {}", size(), i, j, distance);
      return;

    }

    if (j > i) {

      std::swap(i, j);

    }

    packed_matrix_[offset(i, j)] = distance;

  }

  // Tuple returns value, row index, column index in that order.
  [[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> minimum() const {

    if (packed_matrix_.empty()) {

      return { 0.0, 0, 0 };

    }

    auto min_iter = std::ranges::min_element(packed_matrix_);
    auto [row, column] = unpackIndex(static_cast<size_t>(min_iter - packed_matrix_.begin()));
    return { *min_iter, row, column };

  }

  // Tuple returns value, row index, column index in that order.
  [[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> maximum() const {

    if (packed_matrix_.empty()) {

      return { 0.0, 0, 0 };

    }

    auto max_iter = std::ranges::max_element(packed_matrix_);
    auto [row, column] = unpackIndex(static_cast<size_t>(max_iter - packed_matrix_.begin()));
    return { *max_iter, row, column };

  }

  // Rescale elements to the interval [0, 1].
  void normalizeDistance() {

    auto [max, min] = max_min();

    if (max == 0.0) {

      ExecEnv::log().warn("Invalid matrix range; max: {}, min: {}, matrix size: {}", max, min, size());
      return;

    }

    for (size_t row = 0; row < size(); ++row) {
      for (size_t column = 0; column < row; ++column) {

        DistanceType_t raw_distance = getDistance(row, column);
        DistanceType_t adj_distance = raw_distance / max;
        if (not (adj_distance >= 0.0 and adj_distance <= 1.0)) {

          ExecEnv::log().warn("Invalid matrix normalize, raw: {}, range: {}, min: {}, adjusted: {}",
                               raw_distance, max - min, min, adj_distance);

        }
        setDistance(row, column, adj_distance);

      }

    }

  }

private:

  // The .first element is the maximum value, the .second element is the minimum value.
  [[nodiscard]] std::pair<DistanceType_t, DistanceType_t> max_min() const {

    if (packed_matrix_.empty()) {

      return { 0.0, 0.0 };

    }

    auto [min_iter, max_iter] = std::ranges::minmax_element(packed_matrix_);
    return { *max_iter, *min_iter };

  }

  // Packed offset of element (i, j), i > j. Packed order is row-major (matches the reference scan order).
  [[nodiscard]] static size_t offset(size_t i, size_t j) { return packedSize(i) + j; }

  // Number of packed elements below the diagonal of a matrix of the given size.
  [[nodiscard]] static size_t packedSize(size_t matrix_size) { return (matrix_size * (matrix_size - 1)) / 2; }

  // Convert a packed offset back to (row, column) indices, exactly, in integer arithmetic.
  [[nodiscard]] static std::pair<size_t, size_t> unpackIndex(size_t element_offset) {

    size_t row = 1;
    while (packedSize(row + 1) <= element_offset) {

      ++row;

    }
    return { row, element_offset - packedSize(row) };

  }

  size_t matrix_size_{0};
  std::vector<DistanceType_t> packed_matrix_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public class of the UPGMA distance matrix.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


DistanceMatrix::DistanceMatrix() : impl_ptr_(std::make_unique<DistanceMatrixImpl>(1)) {}

DistanceMatrix::DistanceMatrix(size_t matrix_size) : impl_ptr_(std::make_unique<DistanceMatrixImpl>(matrix_size)) {}

DistanceMatrix::DistanceMatrix(DistanceMatrix&& matrix) noexcept {

  impl_ptr_ = std::move(matrix.impl_ptr_);
  matrix.impl_ptr_ = std::make_unique<DistanceMatrixImpl>(1);

}

DistanceMatrix::~DistanceMatrix() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete PIMPL type.

DistanceMatrix& DistanceMatrix::operator=(DistanceMatrix&& matrix) noexcept {

  if (this != &matrix) {

    impl_ptr_ = std::move(matrix.impl_ptr_);
    matrix.impl_ptr_ = std::make_unique<DistanceMatrixImpl>(1);

  }

  return *this;

}


void DistanceMatrix::addMatrix(const DistanceMatrix& add_matrix) { impl().addMatrix(*add_matrix.impl_ptr_); }

void DistanceMatrix::copyMatrix(const DistanceMatrix& copy_matrix) { impl().copyMatrix(*copy_matrix.impl_ptr_); }

[[nodiscard]] size_t DistanceMatrix::size() const { return impl().size(); }

void DistanceMatrix::resize(size_t new_size) { impl().resize(new_size); }

[[nodiscard]] DistanceType_t DistanceMatrix::getDistance(size_t row, size_t column) const { return impl().getDistance(row, column); }

void DistanceMatrix::setDistance(size_t row, size_t column, DistanceType_t distance) { impl().setDistance(row, column, distance); }

// Search efficiency depends on underlying implementation
[[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> DistanceMatrix::minimum() const { return impl().minimum(); }

// Search efficiency depends on underlying implementation
[[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> DistanceMatrix::maximum() const { return impl().maximum(); }

// Rescale elements to the interval [0, 1].
void DistanceMatrix::normalizeDistance() { impl().normalizeDistance(); }


} // Namespace