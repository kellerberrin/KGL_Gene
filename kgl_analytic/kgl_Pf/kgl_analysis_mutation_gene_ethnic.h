//
// Created by kellerberrin on 31/1/21.
//

#ifndef KGL_ANLYSIS_MUTATION_GENE_ETHNIC_H
#define KGL_ANLYSIS_MUTATION_GENE_ETHNIC_H


#include "kgl_genome_genome.h"
#include "kgl_ped_parser.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class GeneEthnicitySex {

public:

  GeneEthnicitySex() = default;
  GeneEthnicitySex(const GeneEthnicitySex &) = default;
  ~GeneEthnicitySex() = default;

  GeneEthnicitySex &operator=(const GeneEthnicitySex &) = default;

  [[nodiscard]] size_t total() const { return total_; }
  [[nodiscard]] const std::map<std::string, size_t> &population() const { return population_; }
  [[nodiscard]] const std::map<std::string, size_t> &superPopulation() const { return super_population_; }
  [[nodiscard]] size_t male() const { return male_; }
  [[nodiscard]] size_t female() const { return female_; }

  bool pedAnalysis(const GenomeId_t& genome_id,
                   size_t count,
                   const std::shared_ptr<const GenomePEDData>& ped_data);

  void updatePopulations(const std::shared_ptr<const GenomePEDData>& ped_data);

  static void writeHeader(const std::shared_ptr<const GenomePEDData>& ped_data,
                          std::ostream& out_file,
                          char output_delimiter);

  void writeOutput(const std::shared_ptr<const GenomePEDData>& ped_data,
                   std::ostream& out_file,
                   char output_delimiter) const;


private:

  // Total number of samples
  size_t total_{0};
  // Population breakdown
  std::map<std::string, size_t> population_;
  // Super Population breakdown
  std::map<std::string, size_t> super_population_;
  // Sex breakdown.
  size_t male_{0};    // Males that have values for this gene.
  size_t female_{0};  // Females that have values for this gene.



};







} // namespacve





#endif //KGL_ANLYSIS_MUTATION_GENE_ETHNIC_H
