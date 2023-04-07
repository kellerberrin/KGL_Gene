//
// Created by kellerberrin on 3/1/21.
//

#ifndef KGL_ANALYSIS_PFEMP_H
#define KGL_ANALYSIS_PFEMP_H


#include "kgl_analysis_virtual.h"
#include "kgl_genome_collection.h"
#include "kgl_pf7_sample_parser.h"
#include "kgl_pf7_fws_parser.h"
#include "kgl_pf7_distance_parser.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class PfEMPAnalysis : public VirtualAnalysis {

  using VariantGeneMap = std::map<std::string, std::pair<std::shared_ptr<const GeneFeature>, std::shared_ptr<ContigDB>>>;

public:

  PfEMPAnalysis() = default;
  ~PfEMPAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return PFEMPANALYSIS_IDENT_; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<PfEMPAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const ActiveParameterList& named_parameters,
                                         const std::shared_ptr<const AnalysisResources>& resource_ptr) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  constexpr static const char PFEMPANALYSIS_IDENT_[]{"PfEMP"};

  std::string ident_work_directory_;
  // Resources
  std::shared_ptr<const Pf7SampleResource> Pf7_sample_ptr_;
  std::shared_ptr<const Pf7FwsResource> Pf7_fws_ptr_;
  std::shared_ptr<const Pf7DistanceResource> Pf7_distance_ptr_;
  constexpr static const char PF3D7_IDENT_[]{"Pf3D7_62"};
  std::shared_ptr<const GenomeReference> genome_3D7_ptr_;     // The Pf7 variant data was aligned on this genome.
  std::shared_ptr<const GenomeCollection> all_reference_genomes_ptr_;
  // Per chromosome VCF files.
  std::shared_ptr<const PopulationDB> all_population_ptr_;
  std::shared_ptr<const PopulationDB> filtered_population_ptr_;  // Only genomes that have passed quality test.
  constexpr static const double MONOCLONAL_FWS_THRESHOLD{0.95};
  std::shared_ptr<const PopulationDB> monoclonal_population_ptr_;  // Only genomes that have passed quality test.

  void performPFEMP1UPGMA();

  // Available Parameters
  constexpr static const char* NEWICK_FILE_ = "NewickFile";
  constexpr static const char* INTRON_FILE_ = "IntronFile";

  // Intron promoter sequences.
  constexpr static const char* I_PROMOTER_ = "TGTATGTG";
  constexpr static const char* I_COMPLEMENT_PROMOTER_ = "ACATACAC";
  constexpr static const char* I_5_PROMOTER_ = "TCATA";

  // Min seq size for UPGMA analysis.
  constexpr static const size_t MIN_SEQUENCE_LENGTH_ = 10;
  // Gene family ident.
  constexpr static const char* PFEMP1_FAMILY_ = "PFEMP1";
  constexpr static const char* RUF6_FAMILY_ = "RUF6";
  constexpr static const char* RIFIN_FAMILY_ = "RIFIN";
  constexpr static const char* STEVOR_FAMILY_ = "STEVOR";
  constexpr static const char* SURFIN_FAMILY_ = "SURFIN";

  // The genes we are interested in.
  VariantGeneMap select_gene_map_;   // Var rifin stevor RUF6
  VariantGeneMap all_gene_map_;     // All genes.

  // File name constants.
  constexpr static const char* NEWICK_{"newick_"};
  constexpr static const char* NEWICK_EXT_{".txt"};
  constexpr static const char* INTRON_{"intron_"};
  constexpr static const char* INTRON_EXT_{".csv"};
  constexpr static const char* VARIANT_COUNT_{"gene_variant_"};
  constexpr static const char* VARIANT_COUNT_EXT_{".csv"};
  constexpr static const char CSV_DELIMITER_ = ',';

  // Return a vector of genes that have a particular text fragment in the description.
  [[nodiscard]] GeneVector getGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr
                                        , const std::string& description_text) const;

  [[nodiscard]] GeneVector getncRNAGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr) const;

  [[nodiscard]] GeneVector proximityGenes(size_t radius,
                                          const std::shared_ptr<const GeneFeature>& target_ptr,
                                          const GeneVector& gene_vector) const;

  // Analyze the introns of Var family genes.
  void varIntron( const GeneVector& gene_vector, const std::string& intron_file) const;

  void geneFamilyUPGMA( const std::shared_ptr<const GenomeReference>& genome_ptr,
                        const GeneVector& gene_vector,
                        const std::string& upgma_file_name,
                        const std::string& family_text) const;

  void checkDistanceMatrix() const;

  [[nodiscard]] VariantGeneMap getSelectGeneMap(const std::shared_ptr<const GenomeReference>& genome_ptr);
  [[nodiscard]] VariantGeneMap getAllGeneMap(const std::shared_ptr<const GenomeReference>& genome_ptr);
  void getGeneVariants(VariantGeneMap& gene_map, const std::shared_ptr<const PopulationDB>& population_ptr);
  void writeGeneResults(const VariantGeneMap& gene_map, const std::string& file_name);

};


} // namespace








#endif //KGL_KGL_ANALYSIS_PFEMP_H
