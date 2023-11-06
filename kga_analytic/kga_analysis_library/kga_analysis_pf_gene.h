//
// Created by kellerberrin on 3/11/23.
//

#ifndef KGA_ANALYSIS_GENE_PF_H
#define KGA_ANALYSIS_GENE_PF_H


#include "kgl_genome_feature.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace


class AnalysisGenePf  {

public:

  AnalysisGenePf() = default;
  ~AnalysisGenePf() = default;

  static void performGeneAnalysis(const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                           const std::string& ident_work_directory_);

private:

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constants for P. Falciparum analysis

  // Available Parameters
  constexpr static const std::string NEWICK_FILE_{"NewickFile"};
  constexpr static const std::string INTRON_FILE_{"IntronFile"};
  constexpr static const std::string NEWICK_{"newick_"};
  constexpr static const std::string NEWICK_EXT_{".txt"};
  constexpr static const std::string INTRON_{"intron_"};
  constexpr static const std::string INTRON_EXT_{".csv"};
  constexpr static const char CSV_DELIMITER_{','};

  // Var intron promoter sequences.
  constexpr static const std::string I_PROMOTER_{"TGTATGTG"};
  constexpr static const std::string I_COMPLEMENT_PROMOTER_{"ACATACAC"};
  constexpr static const std::string I_5_PROMOTER_{"TCATA"};

  // Min seq size for UPGMA analysis.
  constexpr static const size_t MIN_SEQUENCE_LENGTH_{10};
  // Gene family ident.
  constexpr static const std::string PFEMP1_FAMILY_{"PFEMP1"};
  constexpr static const std::string RUF6_FAMILY_{"RUF6"};
  constexpr static const std::string RIFIN_FAMILY_{"RIFIN"};
  constexpr static const std::string STEVOR_FAMILY_{"STEVOR"};
  constexpr static const std::string SURFIN_FAMILY_ {"SURFIN"};
  constexpr static const std::string TRNA_FAMILY_{"TRNA"};
  constexpr static const std::string RIBOSOME_FAMILY_{"RIBOSOM"};

  [[nodiscard]] static GeneVector getGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                          const std::string& desc_uc_text);

  [[nodiscard]] static GeneVector getncRNAGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                                const std::string& desc_uc_text,
                                                size_t max_size = 0);

  [[nodiscard]] static GeneVector proximityGenes( size_t radius,
                                            const std::shared_ptr<const GeneFeature>& target_ptr,
                                            const GeneVector& gene_vector);

  static void varIntron( const GeneVector& gene_vector, const std::string& intron_file_name);

  static void geneFamilyUPGMA( const std::shared_ptr<const GenomeReference>& genome_ptr,
                        const GeneVector& gene_vector,
                        const std::string& newick_file_name,
                        const std::string& family_text);


};



} // Namespace

#endif //KGA_ANALYSIS_GENE_PF_H
