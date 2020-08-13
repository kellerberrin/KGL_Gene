//
// Created by kellerberrin on 3/7/20.
//

#ifndef KGL_ANALYSIS_MUTATION_H
#define KGL_ANALYSIS_MUTATION_H


#include "kgl_analysis_virtual.h"
#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"


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
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  constexpr static const char* REFERENCE_GENOME_ = "GRCh38";
  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";
  constexpr static const char DELIMITER_ = ',';
  std::string output_file_name_;

  constexpr static const double MAX_MAF_ = 0.05; // Set the variant frequency
  constexpr static const ContigOffset_t VARIANT_SPACING_ = 10000;  // Set the variant spacing.

  // The population variant data.
  std::shared_ptr<const GenomeReference> genome_GRCh38_;
  std::shared_ptr<const DiploidPopulation> diploid_population_;
  std::shared_ptr<const UnphasedPopulation> unphased_population_;
  std::shared_ptr<const GenomePEDData> ped_data_;


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


  [[nodiscard]] bool getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters);
  [[nodiscard]] bool hetHomRatio(std::shared_ptr<const DiploidPopulation> population);
  using future_ret_tuple = std::tuple<GenomeId_t, bool, size_t, size_t, double, double>;
  [[nodiscard]] future_ret_tuple processContig(ContigId_t contig_id, std::shared_ptr<const DiploidGenome> genome_ptr);
  [[nodiscard]] std::tuple<bool, double> alleleFrequency_1000Genome(GenomeId_t genome_id, std::shared_ptr<const Variant> variant);
  [[nodiscard]] std::tuple<bool, double> alleleFrequency_SNPdb(GenomeId_t genome_id, std::shared_ptr<const Variant> variant);
  [[nodiscard]] std::tuple<bool, double> alleleFrequency_Gnomad(GenomeId_t genome_id, std::shared_ptr<const Variant> variant_ptr);
  [[nodiscard]] std::tuple<bool, double> processFloatField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name);
  [[nodiscard]] std::tuple<bool, double> processStringField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name);
  [[nodiscard]] std::optional<std::shared_ptr<const Variant>> lookupUnphasedVariant(std::shared_ptr<const Variant> variant_ptr);

// Process Inbreeding coefficient and relatedness using pre-generated allele locus lists.
  bool hetHomRatioLocus(const std::shared_ptr<const DiploidPopulation>& population) const;
// Called by the threadpool for each genome/sample.
  using locus_ret_tuple = std::tuple<GenomeId_t, size_t, size_t, double, double>;
  locus_ret_tuple processLocusContig(const ContigId_t& contig_id,
                                     const std::shared_ptr<const DiploidContig>& contig_ptr,
                                     const std::string& super_population,
                                     const std::shared_ptr<const ContigVariant>& locus_list) const;

// Write the analysis results to a CSV file.
  bool writeResults( const ContigId_t& contig_id,
                     const std::map<GenomeId_t, std::tuple<size_t, size_t, double, double>>& genome_results_map) const;

// Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
  [[nodiscard]] std::shared_ptr<const ContigVariant> getLocusList(const ContigId_t& contig_id,
                                                                  ContigOffset_t spacing,
                                                                  const std::string& super_pop,
                                                                  double min_frequency) const;

// Data check functions (optional).
  void checkPED() const;
  void joinPopulations() const;

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

#endif //KGL_KGL_ANALYSIS_MUTATION_H
