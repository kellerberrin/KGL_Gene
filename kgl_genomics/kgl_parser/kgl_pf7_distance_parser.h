//
// Created by kellerberrin on 3/04/23.
//

#ifndef KGL_PF7_DISTANCE_PARSER_H
#define KGL_PF7_DISTANCE_PARSER_H



#include "kgl_properties_resource.h"
#include "kgl_square_parser.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Gene nomenclature resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using SampleIndexMap = std::map<std::string, size_t>;


class Pf7DistanceResource : public ResourceBase {

public:

  Pf7DistanceResource(std::string identifier,
                      std::shared_ptr<const SampleIndexMap> index_map_ptr,
                      std::shared_ptr<double[]> distance_matrix_ptr,
                      size_t distance_matrix_size_)
      : ResourceBase(ResourceProperties::PF7DISTANCE_RESOURCE_ID_, std::move(identifier)),
        sample_index_map_ptr_(std::move(index_map_ptr)),
        distance_matrix_ptr_(std::move(distance_matrix_ptr)),
        distance_matrix_size_(distance_matrix_size_),
        sample_size_(sample_index_map_ptr_->size()) {}

  ~Pf7DistanceResource() override = default;

  [[nodiscard]] std::shared_ptr<const SampleIndexMap> distanceSamples() const { return sample_index_map_ptr_; }
  [[nodiscard]] double distanceIndex(size_t x_index, size_t y_index) const;
  [[nodiscard]] double distanceIndex(const std::string& sample_x, const std::string& sample_y) const;

private:

  std::shared_ptr<const SampleIndexMap>  sample_index_map_ptr_;
  std::shared_ptr<const double[]> distance_matrix_ptr_;
  const size_t distance_matrix_size_;
  const size_t sample_size_;

  constexpr static const char* NAN_TEXT_{"nan"};

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParsePf7Distance  {

public:

  ParsePf7Distance() {

    // Create valid pointers.
    distance_matrix_ptr_ = std::shared_ptr<double[]>(new double[1]);
    sample_index_map_ptr_ = std::make_shared<SampleIndexMap>();

  }
  ~ParsePf7Distance() = default;

  [[nodiscard]] bool parsePf7Distance(const std::string& matrix_file_name, const std::string& sampleid_file_name);

  [[nodiscard]] std::shared_ptr<SampleIndexMap> getSampleMap() { return sample_index_map_ptr_; }
  [[nodiscard]] std::shared_ptr<double[]> getDistanceMatrix() { return distance_matrix_ptr_; }
  [[nodiscard]] size_t getDistanceMatrixSize() { return distance_matrix_size_; }

private:


  std::shared_ptr<SampleIndexMap>  sample_index_map_ptr_;
  size_t distance_matrix_size_{0};
  std::shared_ptr<double[]> distance_matrix_ptr_;

  constexpr static const char* NAN_TEXT_{"nan"};

};


} // namespace











#endif //KGL_PF7_DISTANCE_PARSER_H
