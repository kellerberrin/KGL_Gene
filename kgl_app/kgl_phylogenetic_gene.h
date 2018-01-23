//
// Created by kellerberrin on 16/01/18.
//

#ifndef KGL_PHYLOGENETIC_GENE_H
#define KGL_PHYLOGENETIC_GENE_H



#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object implements a series of high level functions for the kgl_phylogenetic application.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// pfATP4 drug target ATP4 sodium pump.
#define PFATP4_MINORITY_CONTIG "chr12"
#define PFATP4_MINORITY_GENE "PF3D7_1211900"
#define PFATP4_MINORITY_SEQUENCE "rna_PF3D7_1211900-1"

#define PFATP4_MALAWI_CONTIG "Pf3D7_12_v3"
#define PFATP4_MALAWI_GENE "PF3D7_1211900"
#define PFATP4_MALAWI_SEQUENCE "PF3D7_1211900.1"

// Erythrocyte membrane (Blood Cell Surface) protein (busy Malawi)
#define _BCS_CONTIG "Pf3D7_12_v3"
#define _BCS_GENE "PF3D7_1255200"
#define  BCS_SEQUENCE "PF3D7_1255200.1"


// Mutant rich Rifin (Malawi)
#define RIFIN_3_CONTIG "Pf3D7_01_v3"
#define RIFIN_3_GENE "PF3D7_0101900"
#define RIFIN_3_SEQUENCE "PF3D7_0101900.1"

// Mutant rich Rifin (Malawi)
#define RIFIN_4_CONTIG "Pf3D7_08_v3"
#define RIFIN_4_GENE "PF3D7_0808900"
#define RIFIN_4_SEQUENCE "PF3D7_0808900.1"

// Rifin very similar mutations (Malawi).
#define RIFIN_1_CONTIG "Pf3D7_07_v3"
#define RIFIN_1_GENE "PF3D7_0711700"
#define RIFIN_1_SEQUENCE "PF3D7_0711700.1"

#define RIFIN_2_CONTIG "Pf3D7_07_v3"
#define RIFIN_2_GENE  "PF3D7_0712600"
#define RIFIN_2_SEQUENCE "PF3D7_0712600.1"

// S-Antigen very busy 5-prime and 3-prime regions (Malawi).
#define S_ANTIGEN_CONTIG "Pf3D7_10_v3"
#define S_ANTIGEN_GENE "PF3D7_1035200"
#define S_ANTIGEN_SEQUENCE "PF3D7_1035200.1"

#define ACTIVE_CONTIG RIFIN_4_CONTIG
#define ACTIVE_GENE RIFIN_4_GENE
#define ACTIVE_SEQUENCE RIFIN_4_SEQUENCE




namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class GeneSummary {

public:

  GeneSummary() = default;
  GeneSummary(const GeneSummary&) = default;
  ~GeneSummary() = default;

  GenomeId_t genome;
  std::shared_ptr<DNA5SequenceCoding> prime5_reference;
  std::vector<std::shared_ptr<DNA5SequenceCoding>> prime5_mutant_vec;
  CompareScore_t prime5_score;
  std::shared_ptr<AminoSequence> sequence_ptr;
  std::vector<std::shared_ptr<AminoSequence>> sequence_mutant_vec;
  OffsetVariantMap variant_map;
  CompareScore_t sequence_score;
  std::shared_ptr<DNA5SequenceCoding> prime3_reference;
  std::vector<std::shared_ptr<DNA5SequenceCoding>> prime3_mutant_vec;
  CompareScore_t prime3_score;

};

using GeneSummaryMap = std::map<GenomeId_t, GeneSummary>;


class GeneAnalysis {

public:

  GeneAnalysis() = default;
  virtual ~GeneAnalysis() = default;

  static bool mutateGene(const ContigId_t& contig,
                         const FeatureIdent_t& gene,
                         const FeatureIdent_t& sequence,
                         std::shared_ptr<const PopulationVariant> population_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         const std::string& fasta_filename);

  static bool mutateRegion(const ContigId_t& contig,
                           ContigOffset_t offset,
                           ContigSize_t region_size,
                           std::shared_ptr<const PopulationVariant> population_ptr,
                           std::shared_ptr<const GenomeDatabase> genome_db_ptr);

  static bool mutateGenomeRegion(const GenomeId_t& genome,
                                 const ContigId_t& contig,
                                 const ContigOffset_t offset,
                                 const ContigSize_t region_size,
                                 std::shared_ptr<const PopulationVariant> population_ptr,
                                 std::shared_ptr<const GenomeDatabase> genome_db_ptr);

private:


  constexpr static const ContigSize_t PRIME_REGION_SIZE = 1000;   // Size of the 5 prime and 3 prime regions.

  static bool mutateGenomeGene(const ContigId_t& contig,
                               const FeatureIdent_t& gene,
                               const FeatureIdent_t& sequence,
                               std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                               std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                               GeneSummaryMap& gene_summary_map);


  static bool mutateGenomeRegion(const ContigId_t& contig,
                                 const ContigOffset_t offset,
                                 const ContigSize_t region_size,
                                 std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                 std::shared_ptr<const GenomeDatabase> genome_db_ptr);







};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_PHYLOGENETIC_GENE_H
