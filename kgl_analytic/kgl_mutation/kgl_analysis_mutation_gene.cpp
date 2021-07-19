//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_uniprot_parser.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"
#include "kgl_variant_sort.h"



#include <fstream>
#include <memory_resource>


namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;

void kgl::GenomeMutation::analysisType() {

  std::string member_text;
  switch(gene_membership_) {

    case VariantGeneMembership::BY_SPAN:
      member_text = "gene variant membership by gene span";
      break;

    case VariantGeneMembership::BY_EXON:
      member_text = "gene variant membership by defined exons";
      break;

    default:
    case VariantGeneMembership::BY_ENSEMBL:
      member_text = "gene variant membership by Ensembl gene code (vep field)";
      break;

  }

  ExecEnv::log().info("Gene variant analysis begins with: {}", member_text);

}

// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::genomeAnalysis( const std::vector<std::string>& target_genes,
                                          const std::shared_ptr<const GenomeReference>& genome_ptr,
                                          const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                          const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                                          const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                                          const std::shared_ptr<const EnsemblHGNCResource>& ensembl_nomenclature_ptr)
{

  // Only execute this function once.


  //  kol::PolicyEvidence evidence_policy(kol::GO::getEvidenceType(kol::GO::EvidenceType::EXPERIMENTAL));

  kol::PolicyEvidence evidence_policy;  // all evidence codes.

  const std::vector<kol::GO::EvidenceCode> automated_electronic_evidence{kol::GO::EvidenceCode::IEA};
  evidence_policy.excludeEvidence(automated_electronic_evidence); // exclude IEA

  std::shared_ptr<const kol::TermAnnotation> term_annotation_ptr(std::make_shared<const kol::TermAnnotation>(evidence_policy,
                                                                                                             genome_ptr->geneOntology().getGafRecordVector(),
                                                                                                             kol::AnnotationGeneName::SYMBOLIC_GENE_ID));

  gene_vector_.clear();
  ExecEnv::log().info("Creating Ontology Cache ...");

  OntologyCache ontology_cache(target_genes, term_annotation_ptr, ontology_db_ptr->goGraph());

  ExecEnv::log().info("Ontology Cache Created");

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [offset, gene_ptr] : contig_ptr->getGeneMap()) {

      std::vector<std::string> name_vec;
      gene_ptr->getAttributes().getName(name_vec);
      std::string name;
      std::string gaf_id;

      if (not name_vec.empty()) {

        name = name_vec.front();
        auto [uniprot_id, symbolic_id] = term_annotation_ptr->getGeneIdentifiers(name);
        gaf_id = uniprot_id;

      }

      // If gaf_id_ is empty then try a lookup with the gene id.
      if (gaf_id.empty()) {

        auto [uniprot_id, symbolic_id] = term_annotation_ptr->getGeneIdentifiers(gene_ptr->id());
        gaf_id = uniprot_id;

      }


      auto [hgnc_id, ensembl_id] = getNomenclature( uniprot_nomenclature_ptr, ensembl_nomenclature_ptr, gene_ptr);

      GeneCharacteristic gene_characteristic;
      gene_characteristic.geneDefinition(gene_ptr, genome_ptr->genomeId(), name, hgnc_id, ensembl_id, gaf_id);
      GeneMutation mutation;
      mutation.gene_characteristic = gene_characteristic;
      mutation.clinvar.updateEthnicity().updatePopulations(genome_aux_data);
      mutation.gene_variants.updateLofEthnicity().updatePopulations(genome_aux_data);
      mutation.gene_variants.updateHighEthnicity().updatePopulations(genome_aux_data);
      mutation.gene_variants.updateModerateEthnicity().updatePopulations(genome_aux_data);
      mutation.ontology.processOntologyStats(mutation.gene_characteristic.geneId(), ontology_cache);
      gene_vector_.push_back(mutation);

    } // Gene.

  } // Contig.

  return true;

}


std::tuple<std::string, std::string>
    kgl::GenomeMutation::getNomenclature( const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                                          const std::shared_ptr<const EnsemblHGNCResource>& ensembl_nomenclature_ptr,
                                          const std::shared_ptr<const GeneFeature>& gene_ptr) {


  std::string ensembl_id;

  std::string hgnc_id = gene_ptr->getAttributes().getHGNC();

  // Retrieve Ensembl gene id, if available.
  if (not hgnc_id.empty()) {

    std::vector<std::string> ensembl_vector = uniprot_nomenclature_ptr->HGNCToEnsembl(hgnc_id);

    if (not ensembl_vector.empty()) {

      ensembl_id = ensembl_vector.front();

    }

  }

  ensembl_id = ensembl_nomenclature_ptr->HGNCToEnsembl(hgnc_id);

  return { hgnc_id, ensembl_id};

}


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr,
                                          const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                          const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                          const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  // Count the ethnic samples in the populations.
  ethnic_statistics_.updatePopulations(genome_aux_data);
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    ethnic_statistics_.pedAnalysis(genome_id, 1, genome_aux_data);

  }
  if (not ethnic_statistics_.auditTotals()) {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; problem with ethnic statistics totals");

  }

  ThreadPool thread_pool(ThreadPool::defaultThreads());
  // A vector for futures.
  std::vector<std::future<GeneMutation>> future_vector;
  // Index by Ensembl gene code.
  std::shared_ptr<const EnsemblIndexMap> ensembl_index_map_ptr = VariantSort::ensemblIndex(unphased_population_ptr);
  ExecEnv::log().info("Unphased variants sorted by Ensembl Gene code: {}, Total Unphased Variants: {}",
                      ensembl_index_map_ptr->size(), unphased_population_ptr->variantCount());
  // Queue a thread for each gene.
  for (auto& gene_mutation : gene_vector_) {

    std::future<GeneMutation> future = thread_pool.enqueueTask(&GenomeMutation::geneSpanAnalysis,
                                                               this,
                                                               population_ptr,
                                                               unphased_population_ptr,
                                                               clinvar_population_ptr,
                                                               genome_aux_data,
                                                               ensembl_index_map_ptr,
                                                               gene_mutation);
    future_vector.push_back(std::move(future));

  } // for genes


  std::vector<GeneMutation> gene_vector;
  // Wait on completed threads
  gene_vector.reserve(future_vector.size());
  for (auto& future : future_vector) {

    gene_vector.emplace_back(future.get());

  }

  // todo: This logic is inefficient, the entire gene vector is copied for each VCF file (24 times). Re-design and Re-code.
  gene_vector_ = std::move(gene_vector);

  ExecEnv::log().info("Gene variant Analysis completes, gene count: {}, total gene variants found: {}",
                      gene_vector_.size(), static_cast<size_t>(gene_variant_count_));
  ExecEnv::log().info("Gene variant Analysis statistics, ensembl candidate variants: {}, genome variants checked: {}, genome variants found: {}",
                      static_cast<size_t>(ensembl_variant_count_), static_cast<size_t>(var_checked_count_), static_cast<size_t>(var_found_count_));

  return true;

}


kgl::GeneMutation kgl::GenomeMutation::geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                                         const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                                         const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                                         const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                                         const std::shared_ptr<const EnsemblIndexMap>& ensembl_index_map_ptr,
                                                         GeneMutation gene_mutation) {

  bool contig_data{false};
  EnsemblHashMap ensembl_hash_map;
  ContigOffset_t lower_bound{0};
  ContigOffset_t upper_bound{0};

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

    const ContigId_t& gene_contig_id = gene_mutation.gene_characteristic.contigId();
    auto contig_opt = genome_ptr->getContig(gene_contig_id);
    if (contig_opt) {

      contig_data = true;
      std::shared_ptr<const ContigDB> all_contig_ptr = contig_opt.value();
      if (ensembl_hash_map.empty() and gene_membership_ == VariantGeneMembership::BY_ENSEMBL) {

        getGeneEnsemblHashMap( *ensembl_index_map_ptr,
                               gene_mutation.gene_characteristic,
                               ensembl_hash_map,
                               lower_bound,
                               upper_bound);

      }

      if (not all_contig_ptr->getMap().empty()) {

        // Which method gene membership of variants is determined.
        std::shared_ptr<const ContigDB> gene_variant_ptr;
        switch(gene_membership_) {

          case VariantGeneMembership::BY_SPAN:
            gene_variant_ptr = getGeneSpan(all_contig_ptr, gene_mutation.gene_characteristic);
            break;

          case VariantGeneMembership::BY_EXON:
            gene_variant_ptr = getGeneExon(all_contig_ptr, gene_mutation.gene_characteristic);
            break;

          default:
          case VariantGeneMembership::BY_ENSEMBL: {

            auto contig_ptr = all_contig_ptr->subset(lower_bound, upper_bound);
            gene_variant_ptr = getGeneEnsemblAlt(contig_ptr, ensembl_hash_map, gene_mutation.gene_characteristic);

          }
            break;

        }

        gene_variant_count_ += gene_variant_ptr->variantCount();

        gene_mutation.clinvar.processClinvar( genome_id, gene_contig_id, clinvar_population_ptr, gene_variant_ptr, genome_aux_data);
        gene_mutation.gene_variants.processVariantStats(genome_id, gene_variant_ptr, unphased_population_ptr, genome_aux_data);

      } // contig not empty

    } // if contig

  } // for genome

  // Generate aggregate variant statistics for each gene.
  if (contig_data) {

    if (not gene_mutation.gene_variants.processSummaryStatistics( population_ptr,
                                                                  ethnic_statistics_,
                                                                  gene_mutation.gene_characteristic.geneId())) {

      ExecEnv::log().warn("GenomeMutation::geneSpanAnalysis; problem with processSummaryStatistics, gene: {}",
                          gene_mutation.gene_characteristic.geneId());

    }

  } else {


    gene_mutation.gene_variants.initializeSummaryStatistics( ethnic_statistics_);

  }

  return gene_mutation;

}


// Gets variants over the whole gene span.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneSpan(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneCharacteristic& gene_char) {

  auto subset_contig = contig_ptr->subset(gene_char.geneBegin(), gene_char.geneEnd());

  return subset_contig;

}

// Get variants only occurring within exons for all mRNA sequences.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneExon(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneCharacteristic& gene_char) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(contig_ptr->contigId()));
  std::shared_ptr<const CodingSequenceArray> coding_sequence_array = GeneFeature::getCodingSequences(gene_char.genePtr());

  for (auto const& [sequence_id, sequence_ptr] : coding_sequence_array->getMap()) {

    for (const auto& [cds_id, cds_ptr] : sequence_ptr->getSortedCDS()) {

      gene_contig->merge(contig_ptr->subset(cds_ptr->sequence().begin(), cds_ptr->sequence().end()));

    }

  }

  return gene_contig;

}


// Get variants matching the ensembl.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneEnsembl( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                          const EnsemblIndexMap& ensembl_index_map,
                                                                          const GeneCharacteristic& gene_char) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(gene_char.contigId()));

  if (gene_char.ensemblId().empty()) {

    return  gene_contig;

  }

  auto lower_iterator = ensembl_index_map.lower_bound(gene_char.ensemblId());
  auto const upper_iterator = ensembl_index_map.upper_bound(gene_char.ensemblId());

  while(lower_iterator != upper_iterator) {

    auto const& [ensembl_id, variant_ptr] = *lower_iterator;

    auto offset_opt = contig_ptr->findOffsetArray(variant_ptr->offset());
    if (offset_opt) {

      for (auto const& offset_variant : offset_opt.value()) {

        if (offset_variant->variantHash() == variant_ptr->variantHash()) {

          // Add the unphased population ptr. Note that this loses phasing information (if present)
          // The problem is that the currently (2021) the gnomad 3.1 files have damaged vep fields.
          // So it's either vep or phase information but not both.
          if (not gene_contig->addVariant(variant_ptr)) {

            ExecEnv::log().error( "GenomeMutation::getGeneEnsembl, unable to add variant: {}",
                                  variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

          }

        }

      }

    }

    ++lower_iterator;

  }

  return gene_contig;

}


// Set up Ensembl map.
void kgl::GenomeMutation::getGeneEnsemblHashMap( const EnsemblIndexMap& ensembl_index_map,
                                                 const GeneCharacteristic& gene_char,
                                                 EnsemblHashMap& ensembl_hash_map,
                                                 ContigOffset_t& lower_bound,
                                                 ContigOffset_t& upper_bound) {

  if (gene_char.ensemblId().empty()) {

    lower_bound = 0;
    upper_bound = 0;
    ensembl_hash_map.clear();

  }


  auto lower_iterator = ensembl_index_map.lower_bound(gene_char.ensemblId());
  auto const upper_iterator = ensembl_index_map.upper_bound(gene_char.ensemblId());

  std::shared_ptr<ContigDB> ensembl_contig_ptr(std::make_shared<ContigDB>(gene_char.contigId()));

  while (lower_iterator != upper_iterator) {

    auto const&[ensembl_id, variant_ptr] = *lower_iterator;

    ensembl_hash_map.emplace(variant_ptr->variantHash(), variant_ptr);

    if (not ensembl_contig_ptr->addVariant(variant_ptr)) {

      ExecEnv::log().error( "GenomeMutation::getGeneEnsemblHashMap, unable to add variant: {}",
                            variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

    }

    ++lower_iterator;

  }

  auto [lower, upper] = ensembl_contig_ptr->offsetBounds();

  lower_bound = lower;
  upper_bound = upper;

  ensembl_variant_count_ += ensembl_hash_map.size();

}


// Get variants matching the ensembl.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneEnsemblAlt( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                             const EnsemblHashMap& ensembl_hash_map,
                                                                             const GeneCharacteristic& gene_char) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(gene_char.contigId()));

  for (auto const& [offset, offset_db_ptr] : contig_ptr->getMap()) {

    for (auto const& variant_ptr :  offset_db_ptr->getVariantArray()) {

      ++var_checked_count_;
      auto result = ensembl_hash_map.find(variant_ptr->variantHash());
      if (result != ensembl_hash_map.end()) {

        auto const& [hash, ensembl_variant_ptr] = *result;
        if (not gene_contig->addVariant(ensembl_variant_ptr)) {

          ExecEnv::log().error( "GenomeMutation::getGeneEnsembl, unable to add variant: {}",
                                ensembl_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        }

        ++var_found_count_;

      }

    }

  }

  return gene_contig;

}

