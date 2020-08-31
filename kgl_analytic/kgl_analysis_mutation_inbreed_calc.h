//
// Created by kellerberrin on 27/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_CALC_H
#define KGL_ANALYSIS_MUTATION_INBREED_CALC_H


#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Called by the threadpool for each genome/sample.
struct LocusResults {

  GenomeId_t genome;
  size_t hetero_count{0};
  size_t homo_count{0};
  size_t total_allele_count{0};
  double inbred_allele_sum{0.0};

};

class InbreedingCalculation  {

public:

  InbreedingCalculation() = delete;
  ~InbreedingCalculation() = delete;


  // The inbreeding Algorithms
  [[nodiscard]] static LocusResults processRitlandLocus( const ContigId_t& genome_id,
                                                         const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                         const std::string& super_population_field,
                                                         const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static LocusResults processRitlandPopulation( const GenomeId_t &genome_id,
                                                              const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                              const std::string& super_population_field,
                                                              const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static LocusResults multiLocus1( const ContigId_t& contig_id,
                                                 const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                 const std::string& super_population_field,
                                                 const std::shared_ptr<const ContigVariant>& locus_list);

private:


};



} // namespace















#endif //KGL_ANALYSIS_MUTATION_INBREED_CALC_H
