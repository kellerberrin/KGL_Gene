//
// Created by kellerberrin on 3/1/21.
//

#ifndef KGL_ANALYSIS_PFEMP_H
#define KGL_ANALYSIS_PFEMP_H


#include "kgl_analysis_virtual.h"
#include "kgl_genome_collection.h"



namespace kellerberrin::genome {   //  organization::project level namespace


class PfEMPAnalysis : public VirtualAnalysis {

public:

  PfEMPAnalysis() { reference_genomes_ = std::make_shared<GenomeCollection>(); }
  ~PfEMPAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "PfEMP"; }
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

  std::shared_ptr<GenomeCollection> reference_genomes_;
  std::string ident_work_directory_;

  void performPFEMP1UPGMA();

  // Available Parameters
  constexpr static const char* NEWICK_FILE_ = "NewickFile";
  constexpr static const char* INTRON_FILE_ = "IntronFile";
  constexpr static const char INTRON_DELIMITER_ = ',';

  // Intron promoter sequences.
  constexpr static const char* I_PROMOTER_ = "TGTATGTG";
  constexpr static const char* I_COMPLEMENT_PROMOTER_ = "ACATACAC";
  constexpr static const char* I_5_PROMOTER_ = "TCATA";

  // Min seq size for UPGMA analysis.
  constexpr static const size_t MIN_SEQUENCE_LENGTH_ = 10;
  // Var family ident.
  constexpr static const char* PFEMP1_FAMILY_ = "PFEMP1";

  // File name constants.
  constexpr static const char* NEWICK_{"newick_"};
  constexpr static const char* NEWICK_EXT_{".txt"};
  constexpr static const char* INTRON_{"intron_"};
  constexpr static const char* INTRON_EXT_{".csv"};

  // Return a vector of genes that have a particular text fragment in the description.
  [[nodiscard]] GeneVector getGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr
                                        , const std::string& description_text) const;

  // Analyze the introns of Var family genes.
  void varIntron( const GeneVector& gene_vector, const std::string& intron_file) const;

  void geneFamilyUPGMA( const std::shared_ptr<const GenomeReference>& genome_ptr,
                        const GeneVector& gene_vector,
                        const std::string& upgma_file_name,
                        const std::string& family_text) const;

};


} // namespace








#endif //KGL_KGL_ANALYSIS_PFEMP_H
