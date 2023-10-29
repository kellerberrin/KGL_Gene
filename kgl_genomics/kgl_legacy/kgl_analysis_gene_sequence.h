//
// Created by kellerberrin on 16/01/18.
//

#ifndef KGL_PHYLOGENETIC_GENE_H
#define KGL_PHYLOGENETIC_GENE_H



#include "kgl_sequence_distance.h"
#include "kgl_variant_db_population.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_io_gff_fasta.h"
#include "kgl_genome_collection.h"
#include "kga_analysis_sequence_mutation.h"

#include <memory>
#include <fstream>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object implements a series of high level functions for the kgl_classification application.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// pfATP4 drug target ATP4 sodium pump.

#define PFATP4_CONTIG "Pf3D7_12_v3"
#define PFATP4_GENE "PF3D7_1211900"
#define PFATP4_SEQUENCE "PF3D7_1211900.1"

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




namespace kellerberrin::genome {   //  organization level namespace


class GeneSummary {

public:

  GeneSummary() = default;
  GeneSummary(const GeneSummary&) = default;
  ~GeneSummary() = default;

  GenomeId_t genome;
  std::shared_ptr<DNA5SequenceCoding> prime5_reference;
  std::shared_ptr<DNA5SequenceCoding> prime5_mutant;
  CompareDistance_t prime5_distance;
  std::shared_ptr<AminoSequence> sequence_ptr;
  std::shared_ptr<AminoSequence> sequence_mutant;
  OffsetVariantMap variant_map;
  CompareDistance_t sequence_distance;
  std::shared_ptr<DNA5SequenceCoding> prime3_reference;
  std::shared_ptr<DNA5SequenceCoding> prime3_mutant;
  CompareDistance_t prime3_distance;

};

using GeneSummaryMap = std::map<GenomeId_t, GeneSummary>;


class GenomicSequence {

public:

  GenomicSequence() = delete;
  virtual ~GenomicSequence() = delete;

  [[nodiscard]] static bool mutateGene( const ContigId_t& contig,
                                        const FeatureIdent_t& gene,
                                        const FeatureIdent_t& sequence,
                                        const std::shared_ptr<const PopulationDB>& population_ptr,
                                        const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                        const std::string& fasta_filename);

  [[nodiscard]] static bool mutateGenomeRegion( const GenomeId_t& genome,
                                                const ContigId_t& contig,
                                                const ContigOffset_t offset,
                                                const ContigSize_t region_size,
                                                const std::shared_ptr<const PopulationDB>& population_ptr,
                                                const std::shared_ptr<const GenomeReference>& genome_db_ptr,
                                                const std::string& fasta_file);

  [[nodiscard]] static bool mutateAllRegions( const std::string& file_name,
                                              ContigSize_t region_size,
                                              const std::shared_ptr<const LinearDNASequenceDistance>& dna_distance_metric,
                                              const std::shared_ptr<const PopulationDB>& population_ptr,
                                              const std::shared_ptr<const GenomeReference>& genome_db_ptr);

  [[nodiscard]] static std::string outputRegionHeader(char delimiter);

  [[nodiscard]] static std::string outputGenomeRegion( char delimiter,
                                                       const std::shared_ptr<const LinearDNASequenceDistance>& dna_distance_metric,
                                                       const ContigId_t& contig_id,
                                                       const ContigOffset_t offset,
                                                       const ContigSize_t region_size,
                                                       const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                                       const std::shared_ptr<const GenomeReference>& genome_ref_ptr);

  [[nodiscard]] static bool translateContig( const GenomeId_t& genome_id,
                                             const ContigId_t& gene_id,
                                             const std::shared_ptr<const GenomeCollection>& genome_collection_ptr,
                                             const std::string& fasta_file_name);

  [[nodiscard]] static bool mutateGenomeGene( const ContigId_t& contig_id,
                                              const FeatureIdent_t& gene,
                                              const FeatureIdent_t& sequence,
                                              const std::shared_ptr<const GenomeDB>& genome_variant_ptr,
                                              const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                              GeneSummaryMap& gene_summary_map);

private:


  constexpr static const ContigSize_t PRIME_REGION_SIZE = 1000;   // Default size of the 5 prime and 3 prime regions.



  [[nodiscard]] static bool mutateGenomeRegion( const ContigId_t& contig_id,
                                                const ContigOffset_t offset,
                                                const ContigSize_t region_size,
                                                const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                                const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                                const std::string& fasta_file);



};



}   // end namespace genome



#endif //KGL_PHYLOGENETIC_GENE_H
