//
// Created by kellerberrin on 3/7/20.
//

#ifndef KGL_ANALYSIS_MUTATION_H
#define KGL_ANALYSIS_MUTATION_H


#include "kgl_analysis_virtual.h"
#include "kgl_variant_db_phased.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class MutationAnalysis : public VirtualAnalysis {

public:

  MutationAnalysis() = default;

  ~MutationAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "MUTATION"; }

  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<MutationAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const std::string &work_directory,
                                        const RuntimeParameterMap &named_parameters,
                                        std::shared_ptr<const GenomeCollection> reference_genomes) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const PopulationBase> vcf_data) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  constexpr static const char* REFERENCE_GENOME_ = "GRCh38";
  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";
  constexpr static const char DELIMITER_ = ',';

  std::shared_ptr<const GenomeReference> genome_GRCh38_;
  std::string output_file_name_;
  std::shared_ptr<const DiploidPopulation> filtered_population_;
  std::shared_ptr<const UnphasedPopulation> filtered_joining_population_;

  // Reduce the size of the population to fit into memory.
  ContigOffset_t start_region_ = 0;
  ContigOffset_t end_region_ = 100000000;


  bool getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters);
  bool hetHomRatio(std::shared_ptr<const DiploidPopulation> population);
  using ret_tuple = std::tuple<GenomeId_t, bool, size_t, size_t>;
  ret_tuple processContig(ContigId_t contig_id, std::shared_ptr<const DiploidGenome> genome_ptr);
  void joinPopulations();
};


// Joins a single genome population (Gnomad, Clinvar) to another population.
class JoinSingleGenome {

public:

  JoinSingleGenome(std::shared_ptr<const DiploidPopulation> joined_population,
                   std::shared_ptr<const UnphasedPopulation> joining_population) :
                   joined_population_(joined_population),
                   joining_population_(joining_population) {}
  ~JoinSingleGenome() = default;

  bool joinPopulations();

  // Called by the threadpool per genome.
  std::tuple<std::string, size_t, size_t>
  processGenome(std::shared_ptr<const DiploidGenome> diploid_genome, std::shared_ptr<const GenomeVariant> unphased_genome_ptr);

private:

  // The population to be joined (matched).
  std::shared_ptr<const DiploidPopulation> joined_population_;
  // The population that is matched against the joined population.
  // This is currently assumed to be a single genome population (Gnomad, ExAC, Clinvar, dbSNP).
  std::shared_ptr<const UnphasedPopulation> joining_population_;

};











} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_H
