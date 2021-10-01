//
// Created by kellerberrin on 30/9/21.
//

#ifndef KGL_LITERATURE_DATA_H
#define KGL_LITERATURE_DATA_H

#include <vector>


namespace kellerberrin::genome {   //  organization level namespace



class ImpliedCitation {

public:

  ImpliedCitation() = delete;

  [[nodiscard]] static size_t implied(size_t actual_citations, size_t publication_months);

private:

  // Empirical data for all publications with non-zero citations.
  // 120 month (vector size 120) average citation arrivals in the interval [0.0,1.0].
  // 1.0 indicates all citations over the ten-years.
  // For example, MEAN_CITATION_ARRIVAL[5] = 7.38% is the average proportion of ten years of citations that arrive
  // in the first 6 months after publication.
  static const std::vector<double> MEAN_CITATION_ARRIVAL;

  // The estimated median citation count for research papers.
  const static constexpr size_t MEDIAN_CITATION{30};

};



}

#endif //KGL_LITERATURE_DATA_H
