//
// Created by kellerberrin on 31/1/21.
//

#ifndef KGL_ANLYSIS_MUTATION_GENE_ETHNIC_H
#define KGL_ANLYSIS_MUTATION_GENE_ETHNIC_H


#include "kgl_genome_genome.h"
#include "kgl_hsgenealogy_parser.h"
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

  [[nodiscard]] size_t total() const { return superPopulationTotal(); }
  [[nodiscard]] bool auditTotals() const;
  [[nodiscard]] const std::map<std::string, size_t> &population() const { return population_; }
  [[nodiscard]] const std::map<std::string, size_t> &superPopulation() const { return super_population_; }
  [[nodiscard]] size_t superPopulationCount(const std::string& super_population) const;
  [[nodiscard]] size_t populationCount(const std::string& population) const;
  [[nodiscard]] size_t male() const { return male_; }
  [[nodiscard]] size_t female() const { return female_; }

  // Update population data using genome to lookup the ped record.
  bool genomeAnalysis(const GenomeId_t& genome_id,
                      size_t count,
                      const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  // Generate the map entries for each population
  void updatePopulations(const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  void setDisplay(const std::string& header_prefix, size_t display_flags) {  header_prefix_ = header_prefix; display_flags_ = display_flags; }

  // Write the header for all data, Male, Female, SuperPop, Pop.
  void writeHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                    std::ostream& out_file,
                    char output_delimiter) const;


  void writeOutput(const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                   std::ostream& out_file,
                   char output_delimiter) const;



  constinit const static size_t DISPLAY_SEX_FLAG = 0b1;
  constinit const static size_t DISPLAY_SUPER_POP_FLAG = 0b10;
  constinit const static size_t DISPLAY_POPULATION_FLAG = 0b100;

private:

  // The header prefix.
  std::string header_prefix_{"E_"};
  // The information to output.
  size_t display_flags_{ DISPLAY_SEX_FLAG | DISPLAY_SUPER_POP_FLAG | DISPLAY_POPULATION_FLAG };
  // Population breakdown
  std::map<std::string, size_t> population_;
  // Super Population breakdown
  std::map<std::string, size_t> super_population_;
  // Sex breakdown.
  size_t male_{0};    // Males that have values for this gene.
  size_t female_{0};  // Females that have values for this gene.

  void writeSexHeader( std::ostream& out_file,
                       char output_delimiter) const;

  void writeSuperPopHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                            std::ostream& out_file,
                            char output_delimiter) const;

  void writePopHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                       std::ostream& out_file,
                       char output_delimiter) const;

  void writeSex( std::ostream& out_file,
                 char output_delimiter) const;

  void writeSuperPop( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                      std::ostream& out_file,
                      char output_delimiter) const;

  void writePop( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                 std::ostream& out_file,
                 char output_delimiter) const;

  size_t superPopulationTotal() const;
  size_t populationTotal() const;

};







} // namespacve





#endif //KGL_ANLYSIS_MUTATION_GENE_ETHNIC_H
