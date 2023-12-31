//
// Created by kellerberrin on 30/10/23.
//


#include "kel_exec_env.h"
#include "kgl_distance_matrix_triangular.h"
#include "kgl_distance_matrix.h"


namespace kgl = kellerberrin::genome;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is a fudge around the forward declaration of the PIMPL implementation class.
// There must be a better and more elegant way of achieving this.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DistanceMatrixImplBase = kgl::DistanceMatrixTriangular<kgl::DistanceType_t>;


class kgl::DistanceMatrixImpl : public DistanceMatrixImplBase {

  using DistanceMatrixImplBase::DistanceMatrixImplBase; // Inherit constructors.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public class of the UPGMA calculateDistance matrix.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::DistanceMatrix::DistanceMatrix() : diagonal_impl_ptr_(std::make_unique<DistanceMatrixImpl>(1)) {}

kgl::DistanceMatrix::DistanceMatrix(DistanceMatrix&& matrix) noexcept {

  diagonal_impl_ptr_ = std::move(matrix.diagonal_impl_ptr_);
  matrix.diagonal_impl_ptr_ = std::make_unique<DistanceMatrixImpl>(1);

}

kgl::DistanceMatrix::DistanceMatrix(size_t matrix_size) : diagonal_impl_ptr_(std::make_unique<DistanceMatrixImpl>(matrix_size)) {}

kgl::DistanceMatrix::~DistanceMatrix() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete PIMPL type.


size_t kgl::DistanceMatrix::size() const {

  return diagonal_impl_ptr_->size();

}


void kgl::DistanceMatrix::resize(size_t new_size) {

  diagonal_impl_ptr_->resize(new_size);

}


kgl::DistanceType_t kgl::DistanceMatrix::getDistance(size_t row, size_t column) const {

  return diagonal_impl_ptr_->getDistance( row, column);

}

void kgl::DistanceMatrix::setDistance(size_t row, size_t column, DistanceType_t value) {

  diagonal_impl_ptr_->setDistance( row, column, value);

}

// Search efficiency depends on underlying implementation
std::tuple<kgl::DistanceType_t, size_t, size_t> kgl::DistanceMatrix::minimum() const {

  return diagonal_impl_ptr_->minimum();

}

// Search efficiency depends on underlying implementation
std::tuple<kgl::DistanceType_t, size_t, size_t> kgl::DistanceMatrix::maximum() const {

  return diagonal_impl_ptr_->maximum();

}

// Search efficiency depends on underlying implementation
std::pair<kgl::DistanceType_t, kgl::DistanceType_t> kgl::DistanceMatrix::max_min() const {

  return diagonal_impl_ptr_->max_min();

}

// Rescale elements to the interval [0, 1].
void kgl::DistanceMatrix::normalizeDistance() {

  diagonal_impl_ptr_->normalizeDistance();

}

// Zero all diagonal elements.
void kgl::DistanceMatrix::zeroDiagonal() {

  diagonal_impl_ptr_->zeroDiagonal();

}
