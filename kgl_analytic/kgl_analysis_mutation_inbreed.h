//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_analysis_virtual.h"
#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class InbreedingAnalysis  {

public:

  InbreedingAnalysis( std::shared_ptr<const GenomeReference> genome_GRCh38,
                      std::shared_ptr<const DiploidPopulation> diploid_population,
                      std::shared_ptr<const UnphasedPopulation> unphased_population,
                      std::shared_ptr<const GenomePEDData> ped_data,
                      std::string output_file_name) : genome_GRCh38_(genome_GRCh38),
                                                      diploid_population_(diploid_population),
                                                      unphased_population_(unphased_population),
                                                      ped_data_(ped_data),
                                                      output_file_name_(std::move(output_file_name)) {
    // Clear the data file a
    std::ofstream outfile;
    outfile.open(output_file_name_, std::ofstream::out | std::ofstream::trunc);

  }
  ~InbreedingAnalysis() = default;


  // Data check functions (optional).
  void checkPED() const;
  void joinPopulations() const;

  // Construct a synthetic population and analyze it.
  [[nodiscard]] bool syntheticInbreeding() const;

  // Process Inbreeding coefficient and relatedness using pre-generated allele locus lists.
  bool hetHomRatioLocus() const;

private:

  // The population variant data.
  std::shared_ptr<const GenomeReference> genome_GRCh38_;
  std::shared_ptr<const DiploidPopulation> diploid_population_;
  std::shared_ptr<const UnphasedPopulation> unphased_population_;
  std::shared_ptr<const GenomePEDData> ped_data_;

  constexpr static const char DELIMITER_ = ',';
  std::string output_file_name_;

  constexpr static const double MIN_MAF = 0.01; // Set the minimum MAF
  constexpr static const double MAX_MAF = 0.05; // Set the maximum MAF

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step
  // Synthetic genome contstant
  constexpr static const double SYNTHETIC_GENOME = 1000000; // Used to create the synthetic genome id.


  constexpr static const ContigOffset_t VARIANT_SPACING_ = 1000;  // Set the variant spacing.
  // The Info field identifiers for allele frequency for the 1000 Genome.
  constexpr static const char* GENOME_1000_FREQ_SUFFIX_ = "_AF";
  // The Info field identifiers for allele frequency for db SNP.
  constexpr static const char* SNP_DB_FREQ_ = "FREQ";

  // The Info field identifiers for allele frequency for Gnomad 2.1 and 3.0
  // First is the super populations defined in the 1000 Genomes.
  // Second is the population lookup in the Gnomad Info data.
  constexpr static const std::pair<const char*, const char*> SUPER_POP_AFR_GNOMAD_ {"AFR", "AF_afr"} ;  // African
  constexpr static const std::pair<const char*, const char*> SUPER_POP_AMR_GNOMAD_ {"AMR", "AF_amr"};  // American
  constexpr static const std::pair<const char*, const char*> SUPER_POP_EAS_GNOMAD_ {"EAS", "AF_eas"};  // East Asian
  constexpr static const std::pair<const char*, const char*> SUPER_POP_EUR_GNOMAD_ { "EUR", "AF_nfe"};  // European
  constexpr static const std::pair<const char*, const char*> SUPER_POP_SAS_GNOMAD_ { "SAS", "AF"};  // South Asian
  // Use a super population code to lookup a corresponding AF field.
  std::string lookupSuperPopulationField(const std::string& super_population) const;

  // Synthetic helper functions.
  [[nodiscard]] GenomeId_t generateSyntheticGenomeId(double inbreeding, const std::string& super_population, size_t counter) const;
  [[nodiscard]] std::pair<bool, double> generateInbreeding(const GenomeId_t& genome_id) const;

  [[nodiscard]] bool hetHomRatio(std::shared_ptr<const DiploidPopulation> population) const;
  using future_ret_tuple = std::tuple<GenomeId_t, bool, size_t, size_t, double, double>;
  [[nodiscard]] future_ret_tuple processContig(ContigId_t contig_id, std::shared_ptr<const DiploidGenome> genome_ptr) const;
  [[nodiscard]] std::tuple<bool, double> alleleFrequency_1000Genome(GenomeId_t genome_id, std::shared_ptr<const Variant> variant) const;
  [[nodiscard]] std::tuple<bool, double> alleleFrequency_SNPdb(GenomeId_t genome_id, std::shared_ptr<const Variant> variant) const;
  [[nodiscard]] std::tuple<bool, double> alleleFrequency_Gnomad(GenomeId_t genome_id, std::shared_ptr<const Variant> variant_ptr) const;
  [[nodiscard]] std::tuple<bool, double> processFloatField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) const;
  [[nodiscard]] std::tuple<bool, double> processStringField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) const;
  [[nodiscard]] std::optional<std::shared_ptr<const Variant>> lookupUnphasedVariant(std::shared_ptr<const Variant> variant_ptr) const;

  // Called by the threadpool for each genome/sample.

  struct LocusResults {

    GenomeId_t genome;
    size_t hetero_count{0};
    size_t homo_count{0};
    size_t total_allele_count{0};
    double inbred_allele_sum{0.0};

  };

  LocusResults processRitlandLocus(const ContigId_t& genome_id,
                                   const std::shared_ptr<const DiploidContig>& contig_ptr,
                                   const std::string& super_population_field,
                                   const std::shared_ptr<const ContigVariant>& locus_list) const;

  LocusResults multiLocus1(const ContigId_t& contig_id,
                           const std::shared_ptr<const DiploidContig>& contig_ptr,
                           const std::string& super_population_field,
                           const std::shared_ptr<const ContigVariant>& locus_list) const;

  // Write the analysis results to a CSV file.
  bool writeResults( const ContigId_t& contig_id,
                     const std::map<GenomeId_t, LocusResults>& locus_results) const;

  bool syntheticResults( const ContigId_t& contig_id,
                         const std::map<GenomeId_t, LocusResults>& genome_results_map) const;


  // Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
  // and at a specified frequency for the super population. Used as a template for calculating
  // the inbreeding coefficient and sample relatedness
  [[nodiscard]] std::shared_ptr<const ContigVariant> getLocusList(const ContigId_t& contig_id,
                                                                  ContigOffset_t spacing,
                                                                  const std::string& super_pop,
                                                                  double min_frequency,
                                                                  double max_frequency) const;

  // Create a synthetic population with known inbreeding characteristics
  // Used to test and calibrate the developed inbreeding algorithms.
  std::shared_ptr<const DiploidPopulation> generateSyntheticPopulation( double lower_inbreeding,
                                                                        double upper_inbreeding,
                                                                        double step_inbreeding,
                                                                        const std::string& super_population,
                                                                        const std::shared_ptr<const ContigVariant>& locus_list) const;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Joins a single genome population (Gnomad, Clinvar) to another population.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

class JoinSingleGenome {

public:

  JoinSingleGenome(std::shared_ptr<const DiploidPopulation> joined_population,
                   std::shared_ptr<const UnphasedPopulation> joining_population,
                   std::shared_ptr<const GenomePEDData> ped_data) :
      joined_population_(joined_population),
      joining_population_(joining_population),
      ped_data_(ped_data) {}
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
  // The pedigree data (ethnic background)
  std::shared_ptr<const GenomePEDData> ped_data_;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LocalGenomeJoin {

public :

  explicit LocalGenomeJoin(std::shared_ptr<const GenomeVariant> unphased_genome_ptr) : unphased_genome_ptr_(std::move(unphased_genome_ptr)) {}

  [[nodiscard]] size_t variantsProcessed() const { return variants_processed_; }
  [[nodiscard]] size_t joinedVariantsFound() const { return joined_variants_found_; }

  // Joins a single genome population (Gnomad, Clinvar) to another (generally phased Diploid) population.
  bool lookupJoinedPop(std::shared_ptr<const Variant> variant_ptr);

private:

  size_t variants_processed_{0};
  size_t joined_variants_found_{0};
  std::shared_ptr<const GenomeVariant> unphased_genome_ptr_;


};


} // namespace





#endif //KGL_ANALYSIS_MUTATION_INBREED_H
