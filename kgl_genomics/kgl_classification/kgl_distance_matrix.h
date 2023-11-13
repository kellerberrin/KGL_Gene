//
// Created by kellerberrin on 30/10/23.
//

#ifndef KGL_DISTANCE_MATRIX_H
#define KGL_DISTANCE_MATRIX_H



namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance matrix. Implements the PIMPL pattern to isolate Boost functionality.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using DistanceType_t = double;
class DistanceMatrixImpl;       // Forward declaration of the distance matrix implementation class

class DistanceMatrix {

public:

  DistanceMatrix();
  explicit DistanceMatrix(size_t matrix_size);
  DistanceMatrix(DistanceMatrix&& matrix ) noexcept;
  ~DistanceMatrix();  // Do not use the default destructor, see PIMPL fwd decl below.

  // Tuple returns value, row index, column index in that order.
  [[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> minimum() const;
  [[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> maximum() const;
  // the .first element is the maximum value, the .second element is the minimum value.
  [[nodiscard]] std::pair<DistanceType_t, DistanceType_t> max_min() const;
  // the .first element is the maximum value, the .second element is the minimum value.
  [[nodiscard]] std::pair<DistanceType_t, DistanceType_t> max_min_diagonal() const;


  [[nodiscard]] DistanceType_t getDistance(size_t row, size_t column) const;
  void setDistance(size_t row, size_t column, DistanceType_t distance);

  [[nodiscard]] size_t size() const;
  void resize(size_t new_size);

  // Rescale elements to the interval [0, 1].
  void normalizeDistance();
  // Zero all diagonal elements.
  void zeroDiagonal();

private:

  std::unique_ptr <DistanceMatrixImpl> diagonal_impl_ptr_;    // PIMPL

};


} // Namespace.


#endif //KGL_DISTANCE_MATRIX_H
